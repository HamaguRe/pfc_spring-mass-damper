//! PFCの設計に関わるものをまとめたモジュール

use super::{DVector, DMatrix, Matrix1xX, StateSpace};

pub struct PFC {
    a_m: DMatrix<f64>,  // 内部モデルのシステム行列
    b_m: DMatrix<f64>,  // 内部モデルの入力行列
    x_m: DVector<f64>,  // 内部モデルの状態変数
    k_0: f64,
    nu_x_transpose: Matrix1xX<f64>,
    pub limit: [f64; 2]
}

impl PFC {
    /// * sys: 離散時間状態空間モデル
    /// * n_b: 基底関数の個数
    /// * n_h: 一致点の個数
    /// * t_clrt: 閉ループ応答時間
    /// * limit: 制御入力制約　\[下限, 上限\]
    pub fn new(sys: &StateSpace<f64>, n_b: usize, n_h: usize, t_clrt: f64, limit: [f64; 2]) -> Self {
        let n = sys.a.nrows();
        let (k_0, nu_x) = offline_designer(&sys, n_b, n_h, t_clrt);
        Self {
            a_m: sys.a.clone(),
            b_m: sys.b.clone(),
            x_m: DVector::<f64>::zeros(n),
            k_0: k_0,
            nu_x_transpose: nu_x.transpose(),
            limit: limit
        }
    }

    /// 制御入力を計算して内部モデルを更新する．
    /// 
    /// --- Arguments ---
    /// * r: 目標値
    /// * y: 制御対象出力
    /// 
    /// ---- Return -----
    /// * u: 制御入力
    pub fn update(&mut self, r: f64, y: f64) -> f64 {
        let mut u = self.k_0 * (r - y) + (self.nu_x_transpose.clone() * self.x_m.clone())[0];

        // 入力制約
        if u < self.limit[0] {
            u = self.limit[0]
        } else if u > self.limit[1] {
            u = self.limit[1]
        }

        // 内部モデル更新
        self.x_m = self.a_m.clone() * self.x_m.clone() + self.b_m.clone() * u;
        u
    }
}

/// オフラインでPFCを設計する
/// 
///  ------- Arguments -------
///  * sys: 離散時間状態空間モデル
///  * n_h: 一致点の個数
///  * n_b: 基底関数の個数
///  * t_clrt: 閉ループ応答時間
fn offline_designer(sys: &StateSpace<f64>, n_b: usize, n_h: usize, t_clrt: f64) -> (f64, DVector<f64>) {
    assert!(sys.sample_time > 0.0);
    assert_ne!(n_b, 0);
    assert_ne!(n_h, 0);
    assert!(t_clrt > 0.0);

    // 参照起動の減衰率
    let alpha = (-3.0 * sys.sample_time / t_clrt).exp();

    // nuとnu_xの計算に使用する行列
    let mut tmp0 = DMatrix::<f64>::zeros(n_b, n_h);
    let mut tmp1 = DMatrix::<f64>::zeros(n_b, n_b);
    let mut tmp2 = DMatrix::<f64>::zeros(n_h, sys.a.nrows());
    let mut tmp3 = DVector::<f64>::zeros(n_h);
    for i in 0..n_h {
        // 一致点のサンプル時刻
        let h_time = (t_clrt / (sys.sample_time * (n_h - i) as f64)).floor() as u32;

        let y_b = calc_y_b(&sys, n_b, h_time);
        tmp0.set_column(i, &y_b);
        tmp1 += y_b.clone() * y_b.transpose();
        tmp2.set_row(i, &( sys.c.clone() * sys.a.pow(h_time) - sys.c.clone() ).row(0));
        tmp3[i] = 1.0 - alpha.powi(h_time as i32);
    }
    let nu = tmp0.transpose() * tmp1.try_inverse().unwrap().column(0);
    let k_0 = nu.transpose() * tmp3;
    let nu_x = -tmp2.transpose() * nu;

    (k_0[0], nu_x)
}

/// 一致点における各基底関数に対するモデル出力を
/// まとめたベクトルを計算する．
/// 
///  ------- Arguments -------
///  * sys: 離散時間状態空間モデル
///  * n_b: 基底関数の個数
///  * h_j: 一致点のサンプル時刻
fn calc_y_b(sys: &StateSpace<f64>, n_b: usize, h_j: u32) -> DVector<f64> {    
    let mut y_b = DVector::<f64>::zeros(n_b);
    for l in 0..n_b as u32 {
        let mut y_bl = 0.0;
        let tmp = h_j - 1;
        for q in 0..h_j {
            let a_pow = sys.a.pow(tmp - q);
            let coef = q.pow(l) as f64;
            y_bl += (sys.c.clone() * a_pow * sys.b.clone() * coef)[0];
        }
        y_b[l as usize] = y_bl;
    }
    y_b
}
