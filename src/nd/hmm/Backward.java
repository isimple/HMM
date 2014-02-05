package nd.hmm;

/**
 * 后向算法
 * the same as forward algorithm
 * 
 * @author yang
 */
public class Backward {

	/**
	 * 后向算法内容
	 * 
	 * @param hmm HMM模型
	 * @param o 观察序列
	 * @return 这里返回观察序列的概率
	 */
	public static double standard(HMM hmm, int[] o) {

		double pprob = 0;
		double sum;
		double[][] beta = new double[o.length][hmm.N];

		// 初始化
		for (int i = 0; i < hmm.N; i++) {
			beta[o.length - 1][i] = 1;
		}
		// 迭代计算
		for (int t = o.length - 2; t >= 0; t--) {
			for (int i = 0; i < hmm.N; i++) {
				sum = 0;
				for (int j = 0; j < hmm.N; j++) {
					sum += hmm.A[i][j] * hmm.B[j][o[t + 1]] * beta[t + 1][j];
				}
				beta[t][i] = sum;
			}
		}

		for (int i = 0; i < hmm.N; i++) {
			pprob += beta[0][i];
		}
		return pprob;
	}

	/**
	 * 后向算法估计参数（带比例因子修正）
	 * 
	 * @param hmm HMM模型
	 * @param o 观察序列
	 * @param scale 比例因子
	 * @return 返回概率
	 */

	public static double withScale(HMM hmm, int[] o) {
		double pprob = 0;
		double sum, tmp;
		tmp = 0;
		double[][] beta = new double[o.length][hmm.N];
		double[] scale = new double[o.length];

		// 初始化
		for (int i = 0; i < hmm.N; i++) {
			beta[o.length - 1][i] = 1.0 / scale[o.length - 1];
		}

		// 迭代计算
		for (int t = o.length - 2; t >= 0; t--) {
			for (int i = 0; i < hmm.N; i++) {
				sum = 0;
				for (int j = 0; j < hmm.N; j++) {
					sum += hmm.A[i][j] * hmm.B[j][o[t + 1]] * beta[t + 1][j];
				}
				beta[t][i] = sum / scale[t];
			}
		}

		for (int i = 0; i < hmm.N; i++) {
			tmp += hmm.pi[i] * hmm.B[i][o[0]] * beta[0][i];
		}

		for (int i = 1; i < o.length; i++) {
			pprob += Math.log(scale[i]);
		}
		pprob += Math.log(tmp);
		return pprob;
		
	}

	public static void main(String[] args) {
		double[] pi = new double[] { 0.63, 0.17, 0.2 };
		double[][] A = { { 0.5, 0.375, 0.125 }, { 0.25, 0.125, 0.625 },
				{ 0.25, 0.375, 0.375 } };
		double[][] B = { { 0.6, 0.2, 0.15, 0.05 }, { 0.25, 0.25, 0.25, 0.25 },
				{ 0.05, 0.1, 0.35, 0.5 } };
		HMM hmm = new HMM(3, 4);
		hmm.pi = pi.clone();
		hmm.A = A;
		hmm.B = B;
		int[] T = { 0, 2, 3 };
		//0.026901406250000003
		System.out.println("Forward:" + Forward.standard(hmm, T));
		System.out.println("Backward:" + Backward.standard(hmm, T));
	}

}