package nd.hmm;

public enum HiddenMarkovModelType {

	/**
	 * Specifies a fully connected model,
	 * in which all states are reachable
	 * from all other states.
	 */
	Ergodic,

	/**
	 * Specifies a model in which only forward
	 * state transitions are allowed.
	 */
	Forward,
}
