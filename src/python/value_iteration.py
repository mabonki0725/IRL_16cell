import numpy as np

def trans_mat(env):
    #env.P[s][a][0][1] is next position from s by a action
    #np.eye(1,3,6)[0] = [0,0,0,1,0]
    return (np.array([[np.eye(1, env.nS, env.P[s][a][0][1])[0] for a in range(env.nA)] for s in range(env.nS)]),
            np.array([env.P[s].values()[0][0][2] for s in range(env.nS)]))

def value_iteration(trans_probs, reward, gamma=0.9, epsilon=1e-3):
    """Solving an MDP by value iteration."""
    n_states, n_actions, _ = trans_probs.shape
    U1 = {s: 0 for s in range(n_states)}
    while True:
        U = U1.copy()
        delta = 0
        for s in range(n_states):
            rs = reward[s]	    
            U1[s] = rs + gamma * max([sum([p * U[s1] for s1, p in enumerate(trans_probs[s, a, :])])
                                      for a in range(n_actions)])	
	    """
	    for a in range(n_actions):
		for s1,p in enumerate(trans_probs[s,a,:]):
			vsum = p * U[s1]
			print(s,a,s1,p,U[s1],vsum)
	     """		
            delta = max(delta, abs(U1[s] - U[s]))
	    #print('while',s,U1[s],U[s],delta)	
        if delta < epsilon * (1 - gamma) / gamma:
            return U

def expected_utility(a, s, U, trans_probs):
    """The expected utility of doing a in state s, according to the MDP and U."""
    return sum([p * U[s1] for s1, p in enumerate(trans_probs[s, a, :])])

def best_policy(trans_probs, U):
    """
    Given an MDP and a utility function U, determine the best policy,
    as a mapping from state to action.
    """
    n_states, n_actions, _ = trans_probs.shape
    pi = {}
    for s in range(n_states):
        pi[s] = max(range(n_actions), key=lambda a: expected_utility(a, s, U, trans_probs))
    return pi

if __name__ == '__main__':
    from envs import gridworld
    grid = gridworld.GridworldEnv()
    #print('gird',grid)
    
    fx = np.eye(grid.nS)
    #print(fx)
    
    trans_probs, reward = trans_mat(grid)
    U = value_iteration(trans_probs, reward)
    pi = best_policy(trans_probs, U)
    print('Utility',U)
    print('policy',pi)


    import matplotlib.pyplot as plt

    def to_mat(u, shape):
        dst = np.zeros(shape)
        for k, v in u.iteritems():
            dst[k / shape[1], k % shape[1]] = v
        return dst

    def add_arrow(pi, shape):
        for k, v in pi.iteritems():
            if v == gridworld.UP:
                plt.arrow(k / shape[1], k % shape[1], -0.45, 0, head_width=0.05)
            elif v == gridworld.RIGHT:
                plt.arrow(k / shape[1], k % shape[1], 0, 0.45, head_width=0.05)
            elif v == gridworld.DOWN:
                plt.arrow(k / shape[1], k % shape[1], 0.45, 0, head_width=0.05)
            elif v == gridworld.LEFT:
                plt.arrow(k / shape[1], k % shape[1], 0, -0.45, head_width=0.05)

    plt.matshow(to_mat(U, grid.shape))
    add_arrow(pi, grid.shape)
    plt.show()
