package nd.hmm;

/**
 * HMM对象类
 * 
 * @author yang
 */
public class HMM implements Cloneable {
	/**
	 * 隐藏状态数目
	 */
	public int N;
	/**
	 * 观察状态数目
	 */
	public int M;
	/**
	 * 状态转移矩阵
	 * 一个隐状态到另一个隐状态的概率
	 */
	public double[][] A;

	/**
	 * 混淆矩阵
	 * 一个隐状态到观察状态的概率
	 */
	public double[][] B;
	/**
	 * 初始向量
	 */
	public double[] pi;

	@Override
	public Object clone() {
		HMM hmm = null;
		try {
			hmm = (HMM) super.clone();
			hmm.A = A.clone();
			for (int i = 0; i < A.length; i++) {
				hmm.A[i] = A[i].clone();
			}
			hmm.B = B.clone();
			for (int i = 0; i < B.length; i++) {
				hmm.B[i] = B[i].clone();
			}
			hmm.pi = pi.clone();
		} catch (CloneNotSupportedException e) {
			e.printStackTrace();
		}

		return hmm;
	}

	public HMM() {
		super();
	}

	/**
	 * @param n 隐藏状态数目
	 * @param m 观察状态数目
	 * @param a 状态转移矩阵
	 * @param b 混淆矩阵
	 * @param pi 初始向量
	 */
	public HMM(int n, int m, double[][] a, double[][] b, double[] pi) {
		super();
		N = n;
		M = m;
		A = a.clone();
		for (int i = 0; i < a.length; i++) {
			A[i] = a[i].clone();
		}
		B = b.clone();
		for (int i = 0; i < b.length; i++) {
			B[i] = b[i].clone();
		}
		this.pi = pi.clone();
	}

	/**
	 * 用于参数估计
	 * 
	 * @param n 隐藏状态数目
	 * @param m 观察状态数目
	 */
	public HMM(int n, int m) {
		super();
		N = n;
		M = m;
		A = new double[N][N];
		B = new double[N][M];
		pi = new double[N];
	}

	/**
	 * 用于测试已知模型
	 * @param a 状态转移矩阵
	 * @param b 符号输出矩阵
	 * @param pi 初始向量
	 */
	public HMM(double[][] a, double[][] b, double[] pi) {
		super();
		N = a.length;
		M = b[0].length;
		A = a.clone();
		for (int i = 0; i < a.length; i++) {
			A[i] = a[i].clone();
		}
		B = b;
		for (int i = 0; i < b.length; i++) {
			B[i] = b[i].clone();
		}
		this.pi = pi.clone();
	}


}
