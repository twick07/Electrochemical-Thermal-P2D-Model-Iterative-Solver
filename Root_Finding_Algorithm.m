%Solver With Bisection Method
if solver == 1
    if err(1)*err(3) < 0 

        beta1 = beta1;

        beta2 = beta3;
        
        beta3 = (beta1 +beta2)/2;
     %j1 = beta1*j_ave_neg;
     %j2 = beta2*j_ave_neg;
     %beta3 = (beta1 +beta2)/2;
     %j3 = beta3*j_ave_neg;

    %j_attempts = [j1 j2 j3];
        
        counter_convergence = counter_convergence + 1 ;
    end

    if err(2)*err(3) < 0
        beta1 = beta3;
        beta2 = beta2;
        beta3 = (beta1 +beta2)/2;
        counter_convergence = counter_convergence + 1 ;

     %j1 = beta1*j_ave_neg;
     %j2 = beta2*j_ave_neg;
     %beta3 = (beta1 +beta2)/2;
     %j3 = beta3*j_ave_neg;
    end
end

% Solver with false point method
if solver ==2
    if norm(err(1:3),2) > tol_change
        if err(1)*err(3) < 0 

            beta1 = beta1;

            beta2 = beta3;
        
            beta3 = (beta1 +beta2)/2;

        
        end

        if err(2)*err(3) < 0
            beta1 = beta3;
            beta2 = beta2;
            beta3 = (beta1 +beta2)/2;

        end

    else

    if err(1)*err(3) < 0 

        beta_temp = (beta1*err(2) - beta2*err(1))/(err(2) - err(1)) ;
        beta1 = beta1;

        beta2 = beta3;

        beta3 = beta_temp;
        
        %Mid = (L*calc(U) - U*calc(L))/(calc(U) - calc(L)) ;
        %beta3 = (beta1*err(2) - beta2*err(1))/(err(2) - err(1)) ;

    end
    
    if err(2)*err(3) < 0
        beta_temp = (beta1*err(2) - beta2*err(1))/(err(2) - err(1)) ;
        beta1 = beta3;
        beta2 = beta2;
        beta3 = beta_temp;
    end

    end
    counter_convergence = counter_convergence + 1 ;
end

%Solver for illinois method


if solver == 3
   
    if norm(err(1:3),2) > tol_change
        loop = 1;
        err_prev = err(3);
        if err(1)*err(3) < 0 

            beta1 = beta1;

            beta2 = beta3;
        
            beta3 = (beta1 +beta2)/2;

        
        end

        if err(2)*err(3) < 0
            beta1 = beta3;
            beta2 = beta2;
            beta3 = (beta1 +beta2)/2;

        end
    
    else


    if loop == 1
    if err_prev && err(3) >0
        beta_temp = (0.5*beta1*err(2) - beta2*err(1))/(0.5*err(2) - err(1)) ;

    end

    if err_prev && err(3) < 0
        beta_temp = (beta1*err(2) - 0.5*beta2*err(1))/(err(2) - 0.5*err(1)) ;
    end
    if err(1)*err(3) < 0 
        
        
        beta1 = beta1;

        beta2 = beta3;

        beta3 = beta_temp;
        
        %Mid = (L*calc(U) - U*calc(L))/(calc(U) - calc(L)) ;
        %beta3 = (beta1*err(2) - beta2*err(1))/(err(2) - err(1)) ;

    end
    
    if err(2)*err(3) < 0
        
        beta1 = beta3;
        beta2 = beta2;
        beta3 = beta_temp;
    end
    err_prev = err(3);
    end
    end
    counter_convergence = counter_convergence + 1 ;
end

% solver for Iterate, Truncation and Projection (ITP) method

if solver == 4
    counter_convergence = counter_convergence + 1 ; 
    if norm(err(1:3),2) > tol_change
        if err(1)*err(3) < 0 

            beta1 = beta1;

            beta2 = beta3;
        
            beta3 = (beta1 +beta2)/2;

        
        end

        if err(2)*err(3) < 0
            beta1 = beta3;
            beta2 = beta2;
            beta3 = (beta1 +beta2)/2;

        end

    else
    
    
    
    
    %Calculate All Constants & Interpolate
    U = beta2;
    L = beta1;
    eta_half = log2((U-L)/(2*tol));
    x_bi = (U +L)/2;
    f_U = err(2);
    f_L = err(1);
    x_RF = (U*f_L - L*f_U)/(f_L - f_U) ;
    delta_min = abs(k1*((U-L)^k2));
    delta_max = abs(x_bi - x_RF);
    delta = min([delta_min delta_max]);
    sigma_sol = sign(x_bi - x_RF);

    %Truncation Step
    x_t = x_RF +sigma_sol*delta ; 

    %Projection Step
    row_min = tol*(2^(eta_half +eta_zero - ITP_count)) - (U - L)/2 ;
    row_max = abs(x_t - x_bi) ;
    row = min([row_min,row_max]);
        
    x_ITP = x_bi - sigma_sol*row ;

    ITP_count = ITP_count +1;
    
    %If root is between lower bound and middle
    if err(1)*err(3) < 0

        beta1 = beta1;

        beta2 = x_ITP;

        beta3 = x_ITP;%(beta1 +beta2)/2;


    end

    %If root is between middle and lower boundary
    if err(2)*err(3) < 0
    
        beta1 = x_ITP;

        beta2 = beta2;

        beta3 = x_ITP;%(beta1 +beta2)/2;
    

    end
        
    end  
    







end

%%% Anderson–Björck Method 
if solver == 7
   
    if norm(err(1:3),2) > tol_change
        loop = 1;
        err_prev = err(3);
        if err(1)*err(3) < 0 

            beta1 = beta1;

            beta2 = beta3;
        
            beta3 = (beta1 +beta2)/2;

        
        end

        if err(2)*err(3) < 0
            beta1 = beta3;
            beta2 = beta2;
            beta3 = (beta1 +beta2)/2;

        end
    
    else


    if loop == 1

        %{
            if err_prev && err_next > 0
            var = 1 - calc(Mid)/calc(L) ;

             if var > 0
                m = var;

             else
                
                 m= 0.5;

             end

            Mid = (m*L*calc(U) - U*calc(L))/(m*calc(U) - calc(L)) ;

            err = calc(Mid);
        end

        if err_prev & err_next < 0

             var = 1 - calc(Mid)/calc(U) ;

             if var > 0
                m = var;

             else
                
                 m= 0.5;

             end
            
            Mid = (L*calc(U) - m*U*calc(L))/(calc(U) - m*calc(L)) ;

            err = calc(Mid);

        end
        %}
    if err_prev && err(3) >0

        scale = 1 - err(3)/err(1) ;

        if scale > 0
            scale = scale;
        else
            scale = 0.5;
        end


        beta_temp = (scale*beta1*err(2) - beta2*err(1))/(scale*err(2) - err(1)) ;

    end

    if err_prev && err(3) < 0
        scale = 1 - err(3)/err(2) ;

        if scale > 0
            scale = scale;
        else
            scale = 0.5;
        end

        beta_temp = (beta1*err(2) - scale*beta2*err(1))/(err(2) - scale*err(1)) ;
    end
    if err(1)*err(3) < 0 
        
        
        beta1 = beta1;

        beta2 = beta3;

        beta3 = beta_temp;
        
        %Mid = (L*calc(U) - U*calc(L))/(calc(U) - calc(L)) ;
        %beta3 = (beta1*err(2) - beta2*err(1))/(err(2) - err(1)) ;

    end
    
    if err(2)*err(3) < 0
        
        beta1 = beta3;
        beta2 = beta2;
        beta3 = beta_temp;
    end
    err_prev = err(3);
    end
    end
    counter_convergence = counter_convergence + 1 ;
end

    %mid = U - alpha*calc(U)*(U - L)/(calc(U) - calc(L));
    
    %err = calc(mid);
    
    %L = U;
    %U = mid;