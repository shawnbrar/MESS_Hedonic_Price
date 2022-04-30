function f=f_gmm_sar_he(par,y,x,W1,Q,P1,W)
                [n ~]=size(x);
                b=par(2:end);
                lam=par(1);

                eps = (eye(n)-lam*W1)*y-(x*b);
                P1eps=P1*eps;
                m_t=[eps'*P1eps; Q'*eps];
                f=m_t'*W*m_t;
            end