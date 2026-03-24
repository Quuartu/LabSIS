function K = Cauchy_kernel(u,v,beta)
%CAUCHY_KERNEL computes the kernel matrix corresponding to the input 
%                locations specified by u and v.
%                
%                - u,v are marices whose k-th rows represents the input 
%                  location at the k-th measurement
%                - beta is the hyper-parameter
%                - K is the kernel matrix

K = zeros(size(u,1),size(v,1));
for t=1:size(u,1)
    for s=1:size(v,1)
        K(t,s) = 1/(1+norm(u(t,:)-v(s,:))^2/beta);
    end
end
