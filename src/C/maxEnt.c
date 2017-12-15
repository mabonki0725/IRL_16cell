#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include "mxlib.h"
#include "usrlib.h"

#define GAMMA 0.9
#define EPSILON 0.001
#define ALPHA 0.5

/* 4*4 cell only */
#define CELL16
#define UP     0
#define LEFT   1
#define DOWN   2
#define RIGHT  3

#define N_STATE    16
#define N_ACTION    4
#define N_TRAJS  1000
#define N_FEATURE  16
#define N_T        25 

#define M_REC   4089

typedef struct {
  double prob[N_STATE][N_ACTION][N_STATE];
} TRANS;
typedef struct {
  int ntry;
  int s[N_T];
  int a[N_T];
  int s_next[N_T];
} TRAJ;
typedef struct {
  TRAJ traj[N_TRAJS];
} TRAJS;

/**
  価値関数
**/
int value_iteration(pTrans, reward, U)
TRANS *pTrans;
double *reward;
double *U;  /* N_STATE out */
{
    /* Solving an MDP by value iteration */
    double U1[N_STATE];
    int s,a,s1;
    double delta;
    double Rs;
    double amax,sum;
    double p;

    double limit;
    limit =  EPSILON * (1 - GAMMA) / GAMMA;
    //n_states, n_actions, _ = trans_probs.shape
    //U1 = {s: 0 for s in range(n_states)}
    for(s=0;s<N_STATE;s++) {
      U1[s]=0.0;
      U[s]=0.0;
    }
    while (1) {
      for(s=0;s<N_STATE;s++) U[s] = U1[s];
      delta = 0;
      for(s=0;s<N_STATE;s++) {
        Rs = reward[s];
        amax=-DBL_MAX;
        for(a=0;a<N_ACTION;a++) {
          sum = 0;
          for(s1 =0;s1<N_STATE;s1++) {
            p=pTrans->prob[s][a][s1];
            sum += p * U[s1];
          }
          if(amax < sum) amax = sum;
        }
        //U1[s] = Rs + gamma * max([sum([p * U[s1] for s1, p in enumerate(trans_probs[s, a, :])])
        U1[s] = Rs + GAMMA * amax;
      }      
      delta = -DBL_MAX;
      for(s=0;s<N_STATE;s++) {
        if(delta < fabs(U1[s] - U[s])) delta = fabs(U1[s] - U[s]);
      }
      if(delta < limit) break;    
    }
    return(0);
}
/**
　行動価値関数
**/
double expected_utility(a, s, U, pTrans)
int a;
int s;
double *U;
TRANS *pTrans;
{
    /* The expected utility of doing a in state s, according to the MDP and U. */
    int s1;
    double sum;
    double p;
    //return sum([p * U[s1] for s1, p in enumerate(trans_probs[s, a, :])])
    sum = 0;
    for(s1 =0;s1<N_STATE;s1++) {
      p=pTrans->prob[s][a][s1];
      sum += p * U[s1];
    }
    return(sum);
}
/**
　最適方策
**/
int  best_policy(pTrans, U, pi)
TRANS  *pTrans;
double *U;
int    *pi; /* size state out */
{
    /*""
    Given an MDP and a utility function U, determine the best policy,
    as a mapping from state to action.
    """*/
    int s;
    int a,ma;
    double amax,v;
    double eps;
    
    for(s=0;s<N_STATE;s++) {
       amax=-DBL_MAX;
       ma = 0;
       for(a=0;a<N_ACTION;a++) {
         eps = (double)rand()/(double)RAND_MAX-0.5; /* 摂動を与える */
         v = expected_utility(a,s,U,pTrans);
         if(amax < v*(1+eps/10000.0)) {
           amax=v;
           ma = a;
         }
       }
       pi[s]=ma;
       //pi[s] = max(range(n_actions), key=lambda a: expected_utility(a, s, U, trans_probs))
       //return pi
    }
    return(0);
}
/**
  future
**/
int expected_svf(pTrans, nTrajs, pTrajs, Policy, mu)
TRANS *pTrans;
int nTrajs;
TRAJS *pTrajs;
int *Policy;
double *mu;  /* N_STATE out */
{
    //n_states, n_actions, _ = trans_probs.shape
    //n_t = len(trajs[0])
    int s,t,j,s1;
    double mut[N_STATE][N_T];
    double sum;
    //mu = np.zeros((n_states, n_t))

    for(s=0;s<N_STATE;s++) {
      for(t=0;t<N_T;t++) mut[s][t]=0;
    }
    //for traj in trajs:
    //    mu[traj[0][0], 0] += 1
    //mu[:, 0] = mu[:, 0] / len(trajs)
    for(j=0;j<nTrajs;j++) {
      mut[pTrajs->traj[j].s[0]][0] += 1;
    }
    //for(s=0;s<N_STATE;s++) mut[s][0] /= nTrajs;

    //for t in range(1, n_t):
    //    for s in range(n_states):
    //        mu[s, t] = sum([mu[pre_s, t - 1] * trans_probs[pre_s, policy[pre_s], s] for pre_s in range(n_states)])
    //return np.sum(mu, 1)
#if 0

    for(t=1;t<N_T;t++) {
      for(s=0;s<N_STATE;s++) {
        sum=0;
        for(s1=0;s1<N_STATE;s1++) {
          sum += mut[s1][t-1] * pTrans->prob[s1][Policy[s1]][s];
        }
        mut[s][t]=sum;
      }
    }
#else
    for(j=0;j<nTrajs;j++) {
      for(t=0;t<pTrajs->traj[j].ntry;t++) {
        for(s=0;s<N_STATE;s++) {
          sum=0;
          for(s1=0;s1<N_STATE;s1++) {
            sum += mut[s1][t-1] * pTrans->prob[s1][Policy[s1]][s];
          }
          mut[s][t]=sum;
        }
      }
    }
#endif
    for(s=0;s<N_STATE;s++) {
      for(sum=0,t=0;t<N_T;t++) sum += mut[s][t];
      mu[s] = sum;
    }
    for(s=0;s<N_STATE;s++) mu[s] /= nTrajs;
    return(0);
}
/**
  熟練側
***/
int featureExpert(nTrajs,pTrajs,mu)
int nTrajs;
TRAJS *pTrajs;
double *mu; /* N_STATE out */
{
   int j,s,t;

   for(s=0;s<N_STATE;s++) mu[s]=0;
   for(j=0;j<nTrajs;j++) {
#ifndef CELL16
     s = pTrajs->traj[j].s[0];
     mu[s] += 1.0;
#else
     for(t=0;t<pTrajs->traj[j].ntry;t++) {
       s = pTrajs->traj[j].s[t];
       mu[s] += 1.0;
     }   
#endif
   }
   return(0);
   
}
/**
 IRL
**/       
int max_ent_irl(feature_matrix, pTrans, nTrajs, pTrajs, n_epoch, theta)
double *feature_matrix[];
TRANS *pTrans;
int nTrajs;
TRAJS *pTrajs;
int n_epoch;
double *theta; /* N_FEATURE out */
{
    //n_states, d_states = feature_matrix.shape
    //_, n_actions, _ = trans_probs.shape

    double feature_exp[N_FEATURE];
    int i,L;
    //double theta[N_FEATURE];
    double r[N_STATE];
    double v[N_STATE];
    int    pi[N_STATE];
    double exp_svf[N_STATE];
    double vec[N_FEATURE];
    double grad[N_FEATURE];

    for(i=0;i<N_FEATURE;i++) feature_exp[i]=0.0;
    //feature_exp = np.zeros((d_states))

    //for episode in trajs:
       //for step in episode:
            //feature_exp += feature_matrix[step[0], :]

    featureExpert(nTrajs,pTrajs,feature_exp);
    for(i=0;i<N_FEATURE;i++) {
      feature_exp[i] = feature_exp[i] / nTrajs;
    }


    //theta = np.random.uniform(size=(d_states,))
    for(i=0;i<N_FEATURE;i++) theta[i] = (double)rand()/(double)RAND_MAX;

    for(L=0;L<n_epoch;L++) {
        //r = feature_matrix.dot(theta)
        mxVecR(feature_matrix,N_FEATURE,N_FEATURE,theta,r);
        //v = value_iteration(trans_probs, r, gamma)
        value_iteration(pTrans,r,v);

        best_policy(pTrans, v, pi);

        expected_svf(pTrans, nTrajs, pTrajs, pi, exp_svf);

        mxVecR(feature_matrix,N_FEATURE,N_FEATURE,exp_svf,vec);
        for(i=0;i<N_FEATURE;i++) {
          grad[i] = feature_exp[i] - vec[i];
          theta[i] += ALPHA/100 * grad[i];
        }
    }
    //return feature_matrix.dot(theta)
    return(0);
}
//def feature_matrix(env):
//    return np.eye(env.nS)
/**
　経路の読み込み
**/
int generate_demons(fp, policy, pTrajs)
//int env;
FILE *fp;
double *policy;
TRAJS *pTrajs;
{
  int j,t,k;
  char *p,record[M_REC];
  int gridEnd;
  
  //for _ in range(n_trajs):
  j=0;
#ifdef CELL16
  j=-1;
  while(fgets(record,M_REC,fp)) {
    if(record[0] == '$') continue;
    record[strlen(record)-1]='\0';
    p=strtok(record,", ");
    if(!strcmp(p,"start")) {
      if(j >= 0) { /* 前読込みのPath数 */
        pTrajs->traj[j].ntry=t;
      }
      if(j < N_TRAJS-1) j++; /* 配列超え */
      t=0;
      gridEnd=0;
      continue;
    }
    //if(gridEnd) continue;

    k=0;
    while(p) {
      if(k == 0) pTrajs->traj[j].s[t] = atoi(p);
      if(k == 1) pTrajs->traj[j].a[t] = atoi(p);
      if(k == 2) pTrajs->traj[j].s_next[t] = atoi(p);

      if(k == 2 && (atoi(p) == 0 || atoi(p) == 15)) gridEnd=1; /* 終端判断 */
      k++;
      p=strtok(NULL,", ");
    }
    if(t < N_T-1) t++; /* 配列超え */
  }
#else
  while(fgets(record,M_REC,fp)) {
    if(record[0] == '$') continue;
    record[strlen(record)-1]='\0';
    p=strtok(record,", ");
    t=0;
    k=0;
    while(p) {
      if(k == 0) pTrajs->traj[j].s[t] = atoi(p);
      if(k == 1) pTrajs->traj[j].a[t] = atoi(p);
      if(k == 2) pTrajs->traj[j].s_next[t] = atoi(p);
      k++;
      if(k == 3) {
        t++;
        k = 0;
      }
      if(t >= N_T -1) break; /* 配列超え */

      p=strtok(NULL,", ");
    }
    pTrajs->traj[j].ntry = t;
    j++;
    if(j >= N_TRAJS) break;  /* 配列超え */
  }
#endif
  return(j);
}
/**
  遷移確率
**/
int calcTransProb(nTrajs,pTrajs,pTrans)
int nTrajs;
TRAJS *pTrajs;
TRANS *pTrans;
{
   int j,a,s,t,k,s1;

   for(s=0;s<N_STATE;s++) {
     for(a=0;a<N_ACTION;a++) {
       k=0;
       for(j=0;j<nTrajs;j++) {
         for(t=0;t<pTrajs->traj[j].ntry;t++) {
           if(pTrajs->traj[j].s[t] == s) {
             s1 = pTrajs->traj[j].s_next[t];
             if(pTrajs->traj[j].a[t] == a) pTrans->prob[s][a][s1] += 1;
           }
           k++;
         }
       }
       for(s1=0;s1<N_STATE;s1++) {
         if(k > 0) pTrans->prob[s][a][s1] /= (double)k;
         else      pTrans->prob[s][a][s1] = 0.0;
       }
     }
   }
   return(0);
}
/***
   CELL16 遷移確率
***/
int calcTransCell(pTrans)
TRANS *pTrans;
{
    int a,s,n;
    for(s=0;s<N_STATE;s++) {
      for(a=0;a<N_ACTION;a++) {
        for(n=0;n<N_STATE;n++) pTrans->prob[s][a][n]=0.0;
        if(s == 0 || s==15) pTrans->prob[s][a][s]=1.0;
        else {
          switch(a) {
          case UP   : if(s== 0 || s== 1 || s== 2 || s==3 ) pTrans->prob[s][a][s]=1;
                      else                                 pTrans->prob[s][a][s-4]=1;
                      break;
          case LEFT : if(s== 3 || s== 7 || s==11 || s==15) pTrans->prob[s][a][s]=1;
                      else                                 pTrans->prob[s][a][s+1]=1;
                      break;
          case DOWN : if(s ==12|| s==13 || s==14 || s==15) pTrans->prob[s][a][s]=1;
                      else                                 pTrans->prob[s][a][s+4]=1;
                      break;
          case RIGHT: if(s == 0|| s== 4 || s== 8 || s==12) pTrans->prob[s][a][s]=1;
                      else                                 pTrans->prob[s][a][s-1]=1;
                      break;
          }
        }
      }
    }
    return(0);
}
/**
if __name__ == '__main__':
**/
int main(argc,argv)
int argc;
char *argv[];
{
    //from envs import gridworld
 
    double reward[N_STATE];
    double U[N_STATE];
    int    pi[N_STATE];
    TRANS trans;
    TRAJS trajs;
    double res[N_STATE];
    double **feature_matrix;
    int s;
    int i,j;
    int nTrajs;
    FILE *fp,*fw;

    //grid = gridworld.GridworldEnv()
    //trans_probs, reward = trans_mat(grid)
    if(argc < 3) {
      fprintf(stderr,"USAGE: trajeryFile outptuFile\n");
      return(-9);
    }
    if(!(fp=fopen(argv[1],"r"))) {
      fprintf(stderr,"cannot read trajery file=[%s]\n",argv[1]);
      return(-1);
    }
    if(!(fw=fopen(argv[2],"w"))) {
      fprintf(stderr,"cannot write output file=[%s]\n",argv[2]);
      return(-1);
    }

    feature_matrix=(double **)comMalloc(sizeof(double *)*N_FEATURE);
    for(i=0;i<N_FEATURE;i++) {
      feature_matrix[i]=(double *)comMalloc(sizeof(double)*N_FEATURE);
    }
    
    for(s=0;s<N_STATE;s++) {
      if(s == 0)         reward[s]=0.0;
      else
      if(s == N_STATE-1) reward[s]=0.0;
      else               reward[s]=-1.0;
    }
    //U = value_iteration(trans_probs, reward)
    nTrajs= generate_demons(fp, pi, &trajs);
    fclose(fp);

#ifdef CELL16
    calcTransCell(&trans);  
#else
    calcTransProb(nTrajs,&trajs,&trans);
#endif

    value_iteration(&trans,reward,U);
    best_policy(&trans, U, pi);



    //feature_matrix(grid)
    for(i=0;i<N_FEATURE;i++) {
      for(j=0;j<N_FEATURE;j++) {
        if(i == j) feature_matrix[i][j]=1.0;
        else       feature_matrix[i][j]=0.0;
      }
    }
    max_ent_irl(feature_matrix, &trans, nTrajs, &trajs, 2000, res);

#ifdef CELL16
    for(s=0;s<N_STATE;s++) {
      fprintf(stderr,"reword[%2d]=%lf\n",s,res[s]);
    }
    for(i=0;i<4;i++) {
      for(j=0;j<4;j++) {
        s = i * 4 + j;
        fprintf(fw,"%lf",res[s]);
        if(j < 4-1) fprintf(fw,",");
        else        fprintf(fw,"\n");
      }
    }
#else
    for(s=0;s<N_STATE;s++) {
      fprintf(stderr,"reword[%2d]=%lf\n",s,res[s]);
      fprintf(fw,"%2d,%lf\n",s,res[s]);
    }
#endif

    /* import matplotlib.pyplot as plt
    def to_mat(res, shape):
        dst = np.zeros(shape)
        for i, v in enumerate(res):
            dst[i / shape[1], i % shape[1]] = v
        return dst

    plt.matshow(to_mat(res, grid.shape))
    plt.show()
    **/
    fclose(fw);

    /* free */
    for(i=0;i<N_FEATURE;i++) free(feature_matrix[i]);
    free(feature_matrix);

    return(0);
}