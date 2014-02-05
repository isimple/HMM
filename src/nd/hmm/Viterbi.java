package nd.hmm;

/**
 * Viterbi算法用于计算隐藏状态序列
 * 
 * @author yang
 */
public class Viterbi {

	/**
	 * Viterbi算法
	 * 
	 * @param hmm HMM模型
	 * @param o 观察序列
	 * @param q 最佳隐藏状态序列
	 * @return 最佳路径概率
	 */
	public static double standard(HMM hmm, int[] o, int[] q) {
		double pprob = 0.0;
		//q = new int[o.length];
		double[][] delta = new double[o.length][hmm.N];
		int[][] psi = new int[o.length][hmm.N];

		double maxValue, tmpValue;
		int maxIndex;

		// 初始化
		for (int i = 0; i < hmm.N; i++) {
			delta[0][i] = hmm.pi[i] * hmm.B[i][o[0]];
			psi[0][i] = 0;
		}
		// 迭代计算
		for (int t = 1; t < o.length; t++) {
			for (int j = 0; j < hmm.N; j++) {
				maxIndex = 0;
				maxValue = 0;
				for (int i = 0; i < hmm.N; i++) {
					tmpValue = delta[t - 1][i] * hmm.A[i][j];
					if (tmpValue > maxValue) {
						maxValue = tmpValue;
						maxIndex = i;
					}
				}
				delta[t][j] = maxValue * hmm.B[j][o[t]];
				psi[t][j] = maxIndex;
			}
		}
		// 终止
		q[o.length - 1] = 1;
		for (int i = 0; i < hmm.N; i++) {
			if (delta[o.length - 1][i] > pprob) {
				pprob = delta[o.length - 1][i];
				q[o.length - 1] = i;
			}
		}
		// 求解最佳路径
		for (int t = o.length - 2; t >= 0; t--) {
			q[t] = psi[t + 1][q[t + 1]];
		}
		
		//return pprob;
		return Math.log(pprob);
	}
	
	/**
	 * 取对数方式
	 * @param hmm
	 * @param o
	 * @return
	 */
	public static double withLog(HMM hObj, int[] o,int[] q) {
		HMM hmm = (HMM) hObj.clone();
		double[][] delta = new double[o.length][hmm.N];
		int[][] psi = new int[o.length][hmm.N];
		
		
		double maxValue, tmpValue;
		int maxIndex;
		double pprob = 0;
		double[][] bReplace = new double[hmm.N][o.length];

		// 预处理
		for (int i = 0; i < hmm.N; i++) {
			hmm.pi[i] = Math.log(hmm.pi[i]);
		}
		for (int i = 0; i < hmm.N; i++) {
			for (int j = 0; j < hmm.N; j++) {
				hmm.A[i][j] = Math.log(hmm.A[i][j]);
			}
		}
		for (int i = 0; i < hmm.N; i++) {
			for (int t = 0; t < o.length; t++) {
				bReplace[i][t] = Math.log(hmm.B[i][o[t]]);
			}
		}

		// 初始化
		for (int i = 0; i < hmm.N; i++) {
			//delta[0][i] = hmm.pi[i] * bReplace[i][o[0]];
			delta[0][i] = hmm.pi[i] + bReplace[i][0];
			psi[0][i] = 0;
		}
		// 迭代计算
		for (int t = 1; t < o.length; t++) {
			for (int j = 0; j < hmm.N; j++) {
				maxIndex = 0;
				maxValue = -Double.MAX_VALUE;
				for (int i = 0; i < hmm.N; i++) {
					tmpValue = delta[t - 1][i] + hmm.A[i][j];
					if (tmpValue > maxValue) {
						maxValue = tmpValue;
						maxIndex = i;
					}
				}
				delta[t][j] = maxValue + bReplace[j][t];
				psi[t][j] = maxIndex;
			}
		}
		// 终止
		q[o.length - 1] = 1;
		pprob = -Double.MAX_VALUE;
		for (int i = 0; i < hmm.N; i++) {
			if (delta[o.length - 1][i] > pprob) {
				pprob = delta[o.length - 1][i];
				q[o.length - 1] = i;
			}
		}
		// 求解最佳路径
		for (int t = o.length - 2; t >= 0; t--) {
			q[t] = psi[t + 1][q[t + 1]];
		}

		return pprob;
	}
	
	

	public static void main(String[] args) {
		double[] pi = new double[] { 0.333, 0.333, 0.333 };
		double[][] A = { { 0.333, 0.333, 0.333 }, { 0.333, 0.333, 0.333 },
				{ 0.333, 0.333, 0.333 } };
		double[][] B = { { 0.5, 0.5 }, { 0.75, 0.25 }, { 0.25, 0.75 } };
		HMM hmm = new HMM(A,B,pi);
		int[] T = { 0, 0, 0, 0, 1, 0, 1, 1, 1, 1 };
		//0.026901406250000003
		int[] q1 =  new int[T.length];
		int[] q2 =  new int[T.length];
		System.out.println(Viterbi.standard(hmm, T, q1));
		for (int i = 0; i < q1.length; i++) {
			System.out.print(q1[i] + " ");
		}
		System.out.println("\n**********");
		System.out.println(Viterbi.withLog(hmm, T, q2));
		for (int i = 0; i < q2.length; i++) {
			System.out.print(q2[i] + " ");
		}
		
	}
}
