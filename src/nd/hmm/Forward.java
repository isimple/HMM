package nd.hmm;

/**
 * 前向算法：用于计算观察到序列的概率
 * @author yang
 */
public class Forward {


	/**
	 * 前向算法内容
	 * @param hmm HMM模型
	 * @param o 观察序列
	 * @return 观察序列概率
	 */
	public static double standard(HMM hmm, int[] o) {
		double sum;
		double pprob =0;
		double[][] alpha = new double[hmm.N][hmm.N];
		// 初始化t=1
		for (int i = 0; i < hmm.N; i++) {
			alpha[0][i] = hmm.pi[i] * hmm.B[i][o[0]];
		}
		// 递归
		for (int t = 0; t < o.length - 1; t++) {
			for (int j = 0; j < hmm.N; j++) {
				sum = 0;
				for (int i = 0; i < hmm.N; i++) {
					sum += alpha[t][i] * hmm.A[i][j];
				}
				alpha[t + 1][j] = sum * hmm.B[j][o[t + 1]];

			}
		}
		// 终止
		for (int i = 0; i < hmm.N; i++) {
			pprob += alpha[o.length - 1][i];
		}
		return pprob;
	}
	
	/**
	 * 前向算法估计参数，带比例因子且取对数 解决alpha趋向0，即下溢，除0出现的NaN
	 * @param hmm HMM模型
	 * @param o 观察值序列
	 * @param scale 比例因子
	 * @return
	 */
	public static double withScale(HMM hmm, int[] o) {
		double sum;
		
		double[][] alpha = new double[o.length][hmm.N];
		double[] scale = new double[o.length];
		
		
		// 初始化t=1
		scale[0] = 0;
		for (int i = 0; i < hmm.N; i++) {
			alpha[0][i] = hmm.pi[i] * hmm.B[i][o[0]];
			scale[0] += alpha[0][i];
		}
		for (int i = 0; i < hmm.N; i++) {
			alpha[0][i] /= scale[0];
		}

		// 递归
		for (int t = 0; t < o.length - 1; t++) {
			scale[t + 1] = 0;
			for (int j = 0; j < hmm.N; j++) {
				sum = 0;
				for (int i = 0; i < hmm.N; i++) {
					sum += alpha[t][i] * hmm.A[i][j];
				}
				alpha[t + 1][j] = sum * hmm.B[j][o[t + 1]];
				scale[t + 1] += alpha[t + 1][j];
			}
			for (int i = 0; i < hmm.N; i++) {
				alpha[t + 1][i] /= scale[t + 1];
			}
		}

		double pprob = 0;
		for (int i = 0; i < o.length; i++) {
			pprob  += Math.log10(scale[i]);
		}
		//return pprob;
		return Math.pow(10, pprob);
		
	}
	
	
	
	
	public static void main(String[] args) {
		double[] pi = new double[]{0.63,0.17,0.2};
		double[][] A = {{0.5,0.375,0.125},{0.25,0.125,0.625},{0.25,0.375,0.375}};
		double[][] B = {{0.6,0.2,0.15,0.05},{0.25,0.25,0.25,0.25},{0.05,0.1,0.35,0.5}};
		HMM hmm = new HMM(3,4);
		hmm.pi = pi.clone();
		hmm.A = A;
		hmm.B = B;
		int[] T = {0,2,3};
		//0.026901406250000003
		System.out.println(Forward.standard(hmm, T));
		System.out.println(Forward.withScale(hmm, T));
	}

}
