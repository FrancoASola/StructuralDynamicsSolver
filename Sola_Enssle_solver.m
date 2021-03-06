function [u,ud,udd,fstr,fnd,tan] = Sola_Enssle_solver(mmat,cmat,cnd,kstr,fystr,dtg,ug_dd,dtan)

% Output Variables:
% u,ud,udd are the output vectors of the computed displacement, 
    % velocity and acceleration with size of 3�Nt, where Nt 
    % is the number of discrete time instants  
% f is the output vector of the computed interstory shear from both 
    % columns together. Its size is 3�Nt, where Nt is the number 
    % of discrete time instants
% Input Variables:
% m, c are the 3�3 mass and damping matrices  
% cnd is the constant of the nonlinear dampers (same for both dampers)  
% kstr is the elastic interstory stiffness from both columns together 
    % (same for all stories)  
% fystr is the interstory yield force (same for all stories)  
% dtg is the (scalar) time step of the ground motion,   
% ugdd is the ground acceleration time history in units of g,   
% dtan is the (scalar) time step of the analysis.
err=.000001;                                   % tolerance in norm of R_z vector during NR iterations
jmax=100;                                      % maximum number of NR iterations
ndof=length(mmat);
g=32.2*12;%in/s^2
ug_dd=g*ug_dd;%in/s^2                          % re-define ground acceleration in in/s^2
tg=[0:1:length(ug_dd)-1]*dtg;                  % tg=time vector corresponding to ug_dd (s)
tan=0:dtan:tg(end);                            % tan=vector of time for analysis (s)
ugddan=interp1(tg,ug_dd,tan);                  % ug_ddan = vector of interpolated ground 
                                               % accelerations at analysis time steps
% Initialize marices
u=zeros(ndof,length(tan));
ud=u;
udd=u; udd(:,1)=-ugddan(1);
fstr=u;
Fs=u;
Fd=u;
fnd=u;
Fnd=u;
ktan=zeros(ndof,1);
for i = 1:length(tan)-1
    u_pre=u(:,i);    % u_pre = displacement vector at i (ndofx1)
    ud_pre=ud(:,i);  % ud_pre = velocity vector at i     (ndofx1)
    u_cur=u(:,i);     % u_cur is the ndofx1 vector to become  u(:,i+1) after convergence
    ud_cur=ud(:,i);   % ud_cur is the ndofx1 vector to become ud(:,i+1) after convergence
    dud_du=2/dtan;
%     fprintf('step %g\n',i);
    for j=1:jmax
        %%% Internal column shear force f (ndof x 1)
        fstr(:,i+1)=fstr(:,i)+kstr*[-1  0  0  1  0  0
                               1 -1  0 -1  1  0
                               0  1 -1  0 -1  1]*[u_pre;u_cur];
        for ii=1:ndof
            if abs(fstr(ii,i+1))>fystr
                ktan(ii)=0;
                fstr(ii,i+1)=sign(fstr(ii,i+1))*fystr;
            else
                ktan(ii)=kstr;
            end
        end
        
        ktan_=[ktan(1)+ktan(2)    -ktan(2)          0
               -ktan(2)        ktan(2)+ktan(3)   -ktan(3)
               0                   -ktan(3)       ktan(3)];
           
        %%% Slab Shear Force Fs (ndof x 1)
        Fs(:,i+1)=[fstr(1,i+1)-fstr(2,i+1);
                   fstr(2,i+1)-fstr(3,i+1);
                   fstr(3,i+1)];  %Fs(:,i+1)
               
        %%% Rayleigh (Classical) damping force Fd (ndof x 1)
        Fd(:,i+1)=cmat*ud_cur;
        dFd_du=cmat*dud_du;
        
        %%% Nonlinear damping force Fnd (ndof x 1)
        if abs(ud_cur(1))<=1
            fnd1=cnd*ud_cur(1);
            ctan1=cnd;
        else
            fnd1=cnd*abs(ud_cur(1))^.6*sign(ud_cur(1));
            ctan1=.6*cnd*abs(ud_cur(1))^(-.4);
        end         
        ud_21=ud_cur(2)-ud_cur(1);   % define ud_12 for simplicity and convenience
        if abs(ud_21)<=1
            fnd2=cnd*ud_21;
            ctan2=cnd;
        else
            fnd2=cnd*abs(ud_21)^.6*sign(ud_21);
            ctan2=.6*cnd*abs(ud_21)^(-.4);
        end
        Fnd(:,i+1)=[fnd1-fnd2; fnd2; 0]; %Fnd(:,i+1)
        
        dFnd_du=[ctan1+ctan2        -ctan2         0
                 -ctan2              ctan2         0
                 0                    0            0]*dud_du;
%         R_u_=2*u_i-(2/dtan)*(u_-u_i)+(dtan/2)*m\( ...
%             -Fd(:,i+1)-Fnd(:,i+1)-Fs(:,i+1)-m*ones(ndof,1)*ugddan(i+1) ...
%             -Fd(:,i)-Fnd(:,i)-Fs(:,i)-m*ones(ndof,1)*ugddan(i));

        R_u_cur=(-4/dtan)*mmat*ud_pre+(4/dtan^2)*mmat*(u_cur-u_pre)+ ...
            (Fd(:,i+1)+Fnd(:,i+1)+Fs(:,i+1)+mmat*ones(ndof,1)*ugddan(i+1) ...
            +Fd(:,i)+Fnd(:,i)+Fs(:,i)+mmat*ones(ndof,1)*ugddan(i));
        
        if norm(R_u_cur) < err
            break;
        elseif j>jmax-1
            error(['Newton Raphson Failed to converge after jmax = ' ...
                num2str(jmax) ' iterations.  Fail at i = ' num2str(i)])
        end
%         dR_du=(2/dtan)*eye(3)+(dtan/2)*m\(-ktan_-dFd_du-dFnd_du);
        dR_du=(4/dtan^2)*mmat+(ktan_+dFd_du+dFnd_du);
        u_cur=u_cur-dR_du\R_u_cur;
        ud_cur=-ud_pre+(2/dtan)*(u_cur-u_pre);
    end
fnd(:,i+1)=[fnd1; fnd2; 0];
u(:,i+1)=u_cur;
ud(:,i+1)=ud_cur;
udd(:,i+1)=mmat\(-Fd(:,i+1)-Fnd(:,i+1)-Fs(:,i+1)-mmat*ones(ndof,1)*ugddan(i+1));
end
%end   % function Enssle_Solver