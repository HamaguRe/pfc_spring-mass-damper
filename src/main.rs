//! PFC(Predictive Functional Controller)によるバネ・マス・ダンパ系の位置制御

// 一致点と基底関数の個数によって設計時の行列サイズが変わるのでDMatrixで実装した

use std::fs;
use std::io::{Write, BufWriter};
use nalgebra::{DVector, DMatrix, Matrix1xX};

mod designer;

/// 状態空間モデル
#[derive(Clone)]
pub struct StateSpace<T> {
    a: DMatrix<T>,  // システム行列
    b: DMatrix<T>,  // 入力行列
    c: DMatrix<T>,  // 出力行列
    #[allow(dead_code)]
    d: DMatrix<T>,  // 直達行列
    sample_time: T, // 離散化周期（連続時間なら0にする）
}

impl<T> StateSpace<T> {
    pub fn new(a: DMatrix<T>, b: DMatrix<T>, c: DMatrix<T>, d: DMatrix<T>, dt: T) -> Result<Self, &'static str> {
        let check_n = (a.nrows() == a.ncols()) & (a.nrows() == b.nrows()) & (a.nrows() == c.ncols());
        let check_l = b.ncols() == d.ncols();
        let check_p = c.nrows() == d.nrows();
        if check_n & check_l & check_p {
            Ok(Self {
                a: a,
                b: b,
                c: c,
                d: d,
                sample_time: dt
            })
        } else {
            Err("Matrix size is invalid.")
        }
    }
}

/// ゼロ次ホールドで離散化
/// 
/// * sys: 連続時間状態空間モデル
/// * dt : 離散化周期[s]
fn c2d(sys: StateSpace<f64>, dt: f64) -> StateSpace<f64> {
    assert_eq!(sys.sample_time, 0.0);
    assert!(dt > 0.0);

    let n = sys.a.ncols();
    // 良く見る式
    /*
    let eye = DMatrix::identity(n, n);
    let a_d = (sys.a.clone() * dt).exp();
    // Ax=Bを解く --> A.lu().solve(&B)
    let b_d = sys.a.lu().solve(&((a_d.clone() - eye) * sys.b)).unwrap();
    */
    // こうやってまとめれば一回で計算できる
    let l = sys.b.ncols();
    let mut mat = DMatrix::zeros(n + l, n + l);
    mat.slice_mut((0, 0), (n, n)).copy_from(&(sys.a * dt));
    mat.slice_mut((0, n), (n, l)).copy_from(&(sys.b * dt));
    mat = mat.exp();
    let a_d = mat.slice((0, 0), (n, n)).into();
    let b_d = mat.slice((0, n), (n, l)).into();

    StateSpace::new(a_d, b_d, sys.c, sys.d, dt).unwrap()
}


fn main() {
    // CSVファイルにデータ保存（同一ファイルが存在したら上書き）
    let mut file = BufWriter::new(fs::File::create("result.csv").unwrap());

    // バネ・マス・ダンパ系
    let m = 5.0;  // [kg]
    let c = 5.0;  // [Ns/m]
    let k = 5.0;  // [N/m]
    let plant_c = StateSpace::new(
        DMatrix::from_iterator(2, 2, [
            0.0, 1.0,
            -k/m, -c/m
        ].iter().cloned()).transpose(),  // 正方行列は何故か転置して作られるっぽい
        DMatrix::from_iterator(2, 1, [
            0.0,
            1.0/m
        ].iter().cloned()),
        DMatrix::from_iterator(1, 2, [
            1.0, 0.0
        ].iter().cloned()),
        DMatrix::from_iterator(1, 1, [
            0.0
        ].iter().cloned()),
        0.0
    ).unwrap();

    // 離散化してPFCを設計
    let plant = c2d(plant_c, 0.05);
    let mut pfc = designer::PFC::new(&plant, 2, 3, 0.5, [-5.0, 5.0]);

    let mut x = DVector::zeros(plant.a.nrows());
    for i in 0..=100 {
        let r = if i <= 40 {0.0} else {0.1};
        let y = (plant.c.clone() * x.clone())[0];

        // 制御入力を計算
        let u = pfc.update(r, y);

        // 制御対象の状態を更新
        x = plant.a.clone() * x.clone() + plant.b.clone() * u;

        // データ保存
        file.write(format!(
            "{:.4},{:.4},{:.4},{:.4},{:.4},{:.4}\n",
            plant.sample_time * i as f64, r, y, u, pfc.limit[0], pfc.limit[1]
        ).as_bytes()).unwrap();
    }
}
