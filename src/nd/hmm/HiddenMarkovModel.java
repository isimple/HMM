package nd.hmm;

public class HiddenMarkovModel {

	// Model is defined as M = (A, B, pi)
	/**
	 * Emission probabilities
	 */
	private double[][] B;
	/**
	 * Transition probabilities
	 */
	private double[][] A;
	/**
	 * Initial state probabilities
	 */
	private double[] pi;
	//Size of vocabulary
	private int symbols;
	private int states; // number of states

	public HiddenMarkovModel(int symbols, double[] probabilities)
			throws Exception {
		this(null, probabilities, 0, HiddenMarkovModelType.Ergodic);
		if (symbols <= 0) {
			throw new Exception("Number of symbols should be higher than zero.");
		}

		this.symbols = symbols;
		this.B = new double[states][this.symbols];

		// Initialize B with uniform probabilities
		for (int i = 0; i < this.states; i++)
			for (int j = 0; j < symbols; j++)
				B[i][j] = 1.0 / symbols;
	}

//	public HiddenMarkovModel(double[][] transitions, double[][] emissions,
//			double[] probabilities) {
//		this.symbols = emissions[0].length;
//		this.B = emissions.clone();
//		for (int i = 0; i < emissions.length; i++) {
//			this.B[i] = emissions[i].clone();
//		}
//	}
	public HiddenMarkovModel(double[][] transitions, double[][] emissions,
			double[] probabilities) {
		this.symbols = emissions[0].length;
		this.states = transitions.length;
		this.B = emissions.clone();
		for (int i = 0; i < emissions.length; i++) {
			this.B[i] = emissions[i].clone();
		}
		this.A = transitions.clone();
		for (int i = 0; i < transitions.length; i++) {
			this.A[i] = transitions[i].clone();
		}
		this.pi = probabilities.clone();
		
	}

	public HiddenMarkovModel(int symbols, int states) throws Exception {
		this(symbols, states, HiddenMarkovModelType.Ergodic);
	}

	public HiddenMarkovModel(int symbols, int states, HiddenMarkovModelType type)
			throws Exception {
		this(null, null, states, type);
		if (symbols <= 0) {
			throw new Exception("Number of symbols should be higher than zero.");
		}

		this.symbols = symbols;
		this.B = new double[states][symbols];

		// Initialize B with uniform probabilities
		for (int i = 0; i < states; i++) {
			for (int j = 0; j < symbols; j++) {
				B[i][j] = 1.0 / symbols;
			}
		}
	}

	private HiddenMarkovModel(double[][] transitions, double[] probabilities,
			int states, HiddenMarkovModelType type) throws Exception {
		if (states < 0) {
			throw new Exception("Number of states should be higher than zero.");
		} else if (states == 0) {
			if (probabilities != null) {
				states = probabilities.length;
			} else if (transitions != null) {
				states = transitions.length;
			} else {
				throw new Exception("Number of states could not be determined.");
			}
		}

		int n = states;
		this.states = n;
		if (transitions != null) {
			if (transitions.length != states)
				throw new Exception(
						"Transition matrix should have the same dimensions as the number of model states.");
		} else {

			if (type == HiddenMarkovModelType.Ergodic) {
				// Create A using uniform distribution
				transitions = new double[n][n];
				for (int i = 0; i < n; i++)
					for (int j = 0; j < n; j++)
						transitions[i][j] = 1.0 / n;
			} else {
				// Create A using uniform distribution,
				//   without allowing backward transitions.
				transitions = new double[n][n];
				for (int i = 0; i < n; i++)
					for (int j = i; j < n; j++)
						transitions[i][j] = 1.0 / (n - i);
			}
		}

		this.A = transitions.clone();
		for (int i = 0; i < transitions.length; i++) {
			this.A[i] = transitions[i].clone();
		}
		if (probabilities != null) {
			if (probabilities.length != n)
				throw new Exception(
						"Initial probabilities should have the same length as the number of model states.");
		} else {
			// Create pi as left-to-right
			probabilities = new double[n];
			probabilities[0] = 1.0;
		}

		this.pi = probabilities.clone();

	}

	/**
	 * Calculates the most likely sequence of hidden states
	 * that produced the given observation sequence.
	 * 
	 * @param observations
	 * @param probability new double[1];
	 * @return
	 * @throws Exception
	 */
	public int[] decode(int[] observations, double[] probability)
			throws Exception {
		return decode(observations, false, probability);
	}

	public int[] decode(int[] observations, boolean logarithm,
			double[] probability) throws Exception {
		if (observations == null) {
			throw new Exception("observations");
		}

		if (observations.length == 0) {
			probability[0] = 0.0;
			return new int[0];
		}

		// Viterbi-forward algorithm.
		int T = observations.length;
		int states = this.states;
		int minState;
		double minWeight;
		double weight;

		double[] pi = this.pi;
		double[][] A = this.A;

		int[][] s = new int[states][T];
		double[][] a = new double[states][T];

		// Base
		for (int i = 0; i < states; i++)
			a[i][0] = -Math.log(pi[i]) - Math.log(B[i][observations[0]]);

		// Induction
		for (int t = 1; t < T; t++) {
			int observation = observations[t];

			for (int j = 0; j < states; j++) {
				minState = 0;
				minWeight = a[0][t - 1] - Math.log(A[0][j]);

				for (int i = 1; i < states; i++) {
					weight = a[i][t - 1] - Math.log(A[i][j]);

					if (weight < minWeight) {
						minState = i;
						minWeight = weight;
					}
				}

				a[j][t] = minWeight - Math.log(B[j][observation]);
				s[j][t] = minState;
			}
		}

		// Find minimum value for time T-1
		minState = 0;
		minWeight = a[0][T - 1];

		for (int i = 1; i < states; i++) {
			if (a[i][T - 1] < minWeight) {
				minState = i;
				minWeight = a[i][T - 1];
			}
		}

		// Trackback
		int[] path = new int[T];
		path[T - 1] = minState;

		for (int t = T - 2; t >= 0; t--)
			path[t] = s[path[t + 1]][t + 1];

		// Returns the sequence probability as an out parameter
		probability[0] = logarithm ? -minWeight : Math.exp(-minWeight);

		// Returns the most likely (Viterbi path) for the given sequence
		return path;
	}

	public double evaluate(int[] observations) throws Exception {
		return evaluate(observations, false);
	}

	public double evaluate(int[] observations, boolean logarithm)
			throws Exception {
		if (observations == null)
			throw new Exception("observations");

		if (observations.length == 0)
			return 0.0;

		// Forward algorithm
		double likelihood = 0;
		double[] coefficients = new double[observations.length];

		// Compute forward probabilities
		forward(observations, coefficients);

		for (int i = 0; i < coefficients.length; i++)
			likelihood += Math.log(coefficients[i]);

		// Return the sequence probability
		return logarithm ? likelihood : Math.exp(likelihood);
	}

	public double learn(int[] observations, int iterations, double tolerance)
			throws Exception {
		return learn(new int[][] { observations }, iterations, tolerance);
	}

	public double learn(int[] observations, int iterations) throws Exception {
		return learn(observations, iterations, 0.0);
	}

	public double learn(int[] observations, double tolerance) throws Exception {
		return learn(observations, 0, tolerance);
	}

	public double learn(int[][] observations, double tolerance)
			throws Exception {
		return learn(observations, 0, tolerance);
	}

	public double learn(int[][] observations, int iterations) throws Exception {
		return learn(observations, iterations, 0);
	}

	/**
	 * TODO:the result which input a list of observations
	 * is different from the only observation one 
	 * @param observations
	 * @param iterations
	 * @param tolerance
	 * @return
	 * @throws Exception
	 */
	public double learn(int[][] observations, int iterations, double tolerance)
			throws Exception {
		if (iterations == 0 && tolerance == 0)
			throw new Exception("Iterations and limit cannot be both zero.");

		// Baum-Welch algorithm.

		// The Baum–Welch algorithm is a particular case of a generalized expectation-maximization
		// (GEM) algorithm. It can compute maximum likelihood estimates and posterior mode estimates
		// for the parameters (transition and emission probabilities) of an HMM, when given only
		// emissions as training data.

		// The algorithm has two steps:
		//  - Calculating the forward probability and the backward probability for each HMM state;
		//  - On the basis of this, determining the frequency of the transition-emission pair values
		//    and dividing it by the probability of the entire string. This amounts to calculating
		//    the expected count of the particular transition-emission pair. Each time a particular
		//    transition is found, the value of the quotient of the transition divided by the probability
		//    of the entire string goes up, and this value can then be made the new value of the transition.

		int N = observations.length;
		int currentIteration = 1;
		boolean stop = false;

		double[] pi = this.pi;
		double[][] A = this.A;
		
		// Initialization
		double[][][][] epsilon = new double[N][][][]; // also referred as ksi or psi
		double[][][] gamma = new double[N][][];

		for (int i = 0; i < N; i++) {
			int T = observations[i].length;
			epsilon[i] = new double[T][this.states][this.states];
			gamma[i] = new double[T][this.states];
		}

		// Calculate initial model log-likelihood
		double oldLikelihood = Double.MIN_VALUE;
		double newLikelihood = 0;

		do // Until convergence or max iterations is reached
		{
			// For each sequence in the observations input
			for (int i = 0; i < N; i++) {
				int[] sequence = observations[i];
				int T = sequence.length;
				double[] scaling = new double[T];

				// 1st step - Calculating the forward probability and the
				//            backward probability for each HMM state.
				double[][] fwd = forward(observations[i], scaling); //out 
				double[][] bwd = backward(observations[i], scaling);

				// 2nd step - Determining the frequency of the transition-emission pair values
				//            and dividing it by the probability of the entire string.

				// Calculate gamma values for next computations
				for (int t = 0; t < T; t++) {
					double s = 0;

					for (int k = 0; k < this.states; k++)
						s += gamma[i][t][k] = fwd[t][k] * bwd[t][k];

					if (s != 0) // Scaling
					{
						for (int k = 0; k < states; k++)
							gamma[i][t][k] /= s;
					}
				}

				// Calculate epsilon values for next computations
				for (int t = 0; t < T - 1; t++) {
					double s = 0;

					for (int k = 0; k < states; k++)
						for (int l = 0; l < states; l++)
							s += epsilon[i][t][k][l] = fwd[t][k] * A[k][l]
									* bwd[t + 1][l] * B[l][sequence[t + 1]];

					if (s != 0) // Scaling
					{
						for (int k = 0; k < states; k++)
							for (int l = 0; l < states; l++)
								epsilon[i][t][k][l] /= s;
					}
				}

				// Compute log-likelihood for the given sequence
				for (int t = 0; t < scaling.length; t++)
					newLikelihood += Math.log(scaling[t]);
			}

			// Average the likelihood for all sequences
			newLikelihood /= observations.length;

			// Check if the model has converged or we should stop
			if (checkConvergence(oldLikelihood, newLikelihood,
					currentIteration, iterations, tolerance)) {
				stop = true;
			}

			else {
				// 3. Continue with parameter re-estimation
				currentIteration++;
				oldLikelihood = newLikelihood;
				newLikelihood = 0.0;

				// 3.1 Re-estimation of initial state probabilities 
				for (int k = 0; k < states; k++) {
					double sum = 0;
					for (int i = 0; i < N; i++)
						sum += gamma[i][0][k];
					pi[k] = sum / N;
				}

				// 3.2 Re-estimation of transition probabilities 
				for (int i = 0; i < states; i++) {
					for (int j = 0; j < states; j++) {
						double den = 0, num = 0;

						for (int k = 0; k < N; k++) {
							int T = observations[k].length;

							for (int l = 0; l < T - 1; l++)
								num += epsilon[k][l][i][j];

							for (int l = 0; l < T - 1; l++)
								den += gamma[k][l][i];
						}

						A[i][j] = (den != 0) ? num / den : 0.0;
					}
				}

				// 3.3 Re-estimation of emission probabilities
				for (int i = 0; i < states; i++) {
					for (int j = 0; j < symbols; j++) {
						double den = 0, num = 0;

						for (int k = 0; k < N; k++) {
							int T = observations[k].length;

							for (int l = 0; l < T; l++) {
								if (observations[k][l] == j)
									num += gamma[k][l][i];
							}

							for (int l = 0; l < T; l++)
								den += gamma[k][l][i];
						}

						// avoid locking a parameter in zero.
						B[i][j] = (num == 0) ? 1e-10 : num / den;
					}
				}

			}

		} while (!stop);

		// Returns the model average log-likelihood
		return newLikelihood;
	}

	private double[][] forward(int[] observations, double[] c) {
		int T = observations.length;
		double[] pi = this.pi;
		double[][] A = this.A;
		double[][] fwd = new double[T][states];

		// 1. Initialization
		for (int i = 0; i < states; i++)
			c[0] += fwd[0][i] = pi[i] * B[i][observations[0]];

		if (c[0] != 0) // Scaling
		{
			for (int i = 0; i < states; i++)
				fwd[0][i] = fwd[0][i] / c[0];
		}

		// 2. Induction
		for (int t = 1; t < T; t++) {
			for (int i = 0; i < states; i++) {
				double p = B[i][observations[t]];

				double sum = 0.0;
				for (int j = 0; j < states; j++)
					sum += fwd[t - 1][j] * A[j][i];
				fwd[t][i] = sum * p;

				c[t] += fwd[t][i]; // scaling coefficient
			}

			if (c[t] != 0) // Scaling
			{
				for (int i = 0; i < states; i++)
					fwd[t][i] = fwd[t][i] / c[t];
			}
		}

		return fwd;
	}

	private double[][] backward(int[] observations, double[] c) {
		int T = observations.length;
		//double[] pi = this.pi;
		double[][] A = this.A;

		double[][] bwd = new double[T][states];

		// For backward variables, we use the same scale factors
		//   for each time t as were used for forward variables.

		// 1. Initialization
		for (int i = 0; i < states; i++)
			bwd[T - 1][i] = 1.0 / c[T - 1];

		// 2. Induction
		for (int t = T - 2; t >= 0; t--) {
			for (int i = 0; i < states; i++) {
				double sum = 0;
				for (int j = 0; j < states; j++)
					sum += A[i][j] * B[j][observations[t + 1]] * bwd[t + 1][j];
				bwd[t][i] += sum / c[t];
			}
		}

		return bwd;
	}
	
	public double[][] getB() {
		return B;
	}

	public void setB(double[][] b) {
		B = b;
	}

	public double[][] getA() {
		return A;
	}

	public void setA(double[][] a) {
		A = a;
	}

	public double[] getPi() {
		return pi;
	}

	public void setPi(double[] pi) {
		this.pi = pi;
	}

	public int getSymbols() {
		return symbols;
	}

	public void setSymbols(int symbols) {
		this.symbols = symbols;
	}

	public int getStates() {
		return states;
	}

	public void setStates(int states) {
		this.states = states;
	}

	private static boolean checkConvergence(double oldLikelihood,
			double newLikelihood, int currentIteration, int maxIterations,
			double tolerance) {
		// Update and verify stop criteria
		if (tolerance > 0) {
			// Stopping criteria is likelihood convergence
			if (Math.abs(oldLikelihood - newLikelihood) <= tolerance)
				return true;

			if (maxIterations > 0) {
				// Maximum iterations should also be respected
				if (currentIteration >= maxIterations)
					return true;
			}
		} else {
			// Stopping criteria is number of iterations
			if (currentIteration == maxIterations)
				return true;
		}

		// Check if we have reached an invalid state
		if (Double.isNaN(newLikelihood) || Double.isInfinite(newLikelihood)) {
			return true;
		}

		return false;
	}

	public static void main(String[] args) throws Exception {
		int[][] sequences = new int[][] { new int[] { 0, 1, 1, 1, 1, 1, 1 },
				new int[] { 0, 1, 1, 1 }, new int[] { 0, 1, 1, 1, 1 },
				new int[] { 0, 1, }, new int[] { 0, 1, 1 }, };
		// Creates a new Hidden Markov Model with 2 states for
		//  an output alphabet of two characters (zero and one)
		HiddenMarkovModel hmm = new HiddenMarkovModel(2, 2);

		// Try to fit the model to the data until the difference in
		//  the average likelihood changes only by as little as 0.01
//		for (int i = 0; i < sequences.length; i++) {
//			hmm.learn(sequences[i], 0.01);
//		}
		hmm.learn(sequences, 0.01);
		// Calculate the probability that the given
		//  sequences originated from the model
		double l1 = hmm.evaluate(new int[] { 0, 1 }); // l1 = 0.9999
		double l2 = hmm.evaluate(new int[] { 0, 1, 1, 1 }); // l2 = 0.9999

		double l3 = hmm.evaluate(new int[] { 1, 1 }); // l3 = 0.0000
		double l4 = hmm.evaluate(new int[] { 1, 0, 0, 0 }); // l4 = 0.0000
		System.out.println(l1 + "\t" + l2 + "\t" + l3 + "\t" + l4);
		
		double[] pi = new double[]{0.63,0.17,0.2};
		double[][] A = {{0.5,0.375,0.125},{0.25,0.125,0.625},{0.25,0.375,0.375}};
		double[][] B = {{0.6,0.2,0.15,0.05},{0.25,0.25,0.25,0.25},{0.05,0.1,0.35,0.5}};
		HiddenMarkovModel hmm2 = new HiddenMarkovModel(4,3);
		hmm2.setPi(pi);
		hmm2.setA(A);
		hmm2.setB(B);
		int[] T = {0,2,3};
		//0.026901406250000003
		System.out.println(hmm2.evaluate(T));
		
		double[] pi2 = new double[] { 0.333, 0.333, 0.333 };
		double[][] A2 = { { 0.333, 0.333, 0.333 }, { 0.333, 0.333, 0.333 },
				{ 0.333, 0.333, 0.333 } };
		double[][] B2 = { { 0.5, 0.5 }, { 0.75, 0.25 }, { 0.25, 0.75 } };
		HiddenMarkovModel hmm3 = new HiddenMarkovModel(A2 ,B2,pi2);
		int[] T2 = { 0, 0, 0, 0, 1, 0, 1, 1, 1, 1 };
		//0.026901406250000003
		double[] q1 =  new double[1];
		int[] kk = hmm3.decode(T2, true,q1);
		System.out.println(q1[0]);
		for (int i = 0; i < kk.length; i++) {
			System.out.print(kk[i] + " ");
		}
		
		
	}
}
