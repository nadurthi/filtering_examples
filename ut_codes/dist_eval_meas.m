function f=dist_eval_meas(prior_GMM,z,h,g,R,s,N)
postGMM=post_GMM(prior_GMM,z,h,g,R,s);
f=dist_GM(postGMM,prior_GMM)^N;
end