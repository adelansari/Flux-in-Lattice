%%   Nuclear Reactor Theory                     %%
%*   Matrix assignment (Flux in lattice)        *%
%*   Ahmed Mustafa Ishag       -- U00038807     *%
%*   Adel Ali Ansari           -- U00036658     *%
%************************************************%
clear all
clc

%%  The input code

disp('****Nuclear Reactor Theory');
disp(sprintf('****Matrix assignment (Flux in lattice)\n'));
disp(sprintf('**** Student names and ID numbers:'));
disp(sprintf('\t\t\t\tAhmed Mustafa Ishag       -- U00038807\n'));
disp(sprintf('\t\t\t\tAdel Ali Ansari           -- U00036658\n'));
disp(sprintf('\n\t\t\t\tInstructor:\tDr.Iyad Al Qasir\n'));


%   The number of meshes
slab=input('\n\nSpecify the length of the region: ');
%   The number of meshes
m=input('\n\nSpecify the number of divisions: ');

%divisions
delta=(slab/m);

sigmaA =input('\n\n Enter the macroscopic absorbtion cross section: ');

Diffusion_coefficient =input('\n\n Enter the Diffusion coefficient: ');

a_material= (-Diffusion_coefficient)/(delta^2);
b_material= sigmaA + ((2*Diffusion_coefficient)/(delta^2));

%   Entering the flux values
R=eye(m);

R= b_material*R ;


for i=1:m-1
        R(i,i+1)= a_material;
        
        R(i+1,i)=R(i,i+1);
end

%   displaying the flux coefficients in Matrix form
disp(sprintf('\nFlux Coefficients=\n'));
R

%  Entering the source values
choice=input('\n\n If you would like to manually input the source \n at each division, press 1. \n Alternatively, press 2 to input the source as a function: ');
if choice==1
for i=1:m
    V(i,1)=input(sprintf('\nEnter the source at the division (%d)= ',i));
end

else if choice==2
    rad=1:m;
    V(rad,1)= input('\n\n Specify the source-function in terms of the variable "rad": ') ;
    end
end
disp(sprintf('\nSource=\n'));
S=V

%%exit or continue
conex=input(sprintf('\n\nPress 1 and enter to continue. \nPress 0 and Enter to stop the program. \n1 or 0? '));
while conex ~=0
    
%   Asking the user to choose the method to use for solving the matrix.

disp(sprintf('\n\nPlease specify the method to solve the non-homogenous linear system:'));
disp(sprintf('[1]   Triangular Factorization Method'));
disp(sprintf('[2]   Gaussian Eliminatoin Method'));
disp(sprintf('[3]   Jacobi Iterative method'));
disp(sprintf('[4]   Gauss-Seidel Iterative Method'));
disp(sprintf('[5]   Successive Relaxation Method'));

matsol=input('Enter the number (without brackets): ');


if matsol==1
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%     Triangular Factorization Method       %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %   Test

    C=abs(R);
    solution=1;
    for j=1:m
        b = max(C(j:m,j));
        p = find(C(:,j)==b);
        if C(p,j)==0
            disp(sprintf('\n\t**Error!!! The maximum element in column (%d) is zero.\n',j));
            disp(sprintf('\t  The matrix is singular. There''s no exact solution.\n'));
            solution=0;
            break
        end
    end

    if solution==1
        disp(sprintf('\n\n**There''s an exact solution since matrix R is a nonsingular matrix.\n'));
    end


    %   Initialize L to the identity matrix:
    L = eye(m);

    %   Initialize U to the R matrix:
    U = R;

    for j=1:m-1
        for i=j+1:m
            L(i,j)= U(i,j)/U(j,j);
            U(i,:)=U(i,:)-U(i,j)/U(j,j)*U(j,:);
        end
    end


    %   Finding the matrix Y using forward substitution
    Y= zeros(m,1);   %init vector to zeros
    Y(1)= V(1);
    for i=2:m
        Y(i)=V(i)-L(i,:)*Y;
    end

    %   Finding Flux values, I, using backward substitution
    I= zeros(m,1);
    I(m)=Y(m)/U(m,m);
    for i=m-1:-1:1
        I(i)=(Y(i)-U(i,:)*I)/U(i,i);
    end
    
    %   displaying the Lower Triangular matrix (L)
    disp(sprintf('\nThe Lower Triangular matrix (L):\n'));
    L
    
    %   Displaying the Upper Triangular matrix (U)
    disp(sprintf('\nThe Upper Triangular matrix (U):\n'));
    U
    
    %   Printing the flux values, I
    disp(sprintf('\nThe values of the flux are phi='))
    flux=I
    
    



elseif matsol==2   
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%     Gaussian Eliminatoin Method      %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %  The Algorithm
    C=[R,V]; % The augmented matrix

    solution=1;

    for j=1:m
        p=j; % Finding the row of the maximum element in the column
        for i=j:m-1
            if abs(C(i+1,j))> abs(C(i,j))
                p=i+1;
            end
        end

        %  Checking if there is solution
        if C(p,j)==0
            disp(sprintf('\n\t**Error!The maximum element in column (%d) is zero.\n\t  There''s no exact solution.\n\t  End of the programme.',j));
            solution=0;
            break
        end

        %  Interchanging the rows
        if p>j
            Z=C(p,:);
            C(p,:)=C(j,:);
            C(j,:)=Z;
        end

        for i=j+1:m
            C(i,:)=C(i,:)-(C(i,j)/C(j,j))*C(j,:);
        end

    end

    %  If there is exact solution
    if solution==1
        %  Partition matrix C to matrix D & vector E
        for j=1:m
            for i=1:m
                D(j,i)=C(j,i);
            end
            E(j,1)=C(j,m+1);
        end

        %  Solving flux vector in backward
        I_guess(m,1)=E(j,1)/D(j,j);
        for j=m-1:-1:1
            S=0;
            for i=j+1:m
                S=S+(D(j,i)*I_guess(i,1));
            end
            I_guess(j,1)=(1/D(j,j))*(E(j,1)-S);
        end

        %  Sending the flux obtained by the Gaussian Elimination
        disp(sprintf('\n\nThe flux obtained from the Gaussian Elimination phi\n'));

        flux=I_guess
    end
    
elseif matsol==3
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%     Jacobi Iterative method     %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %   Entering the value of the tolerance
        tol=input(sprintf('\n\nEnter the tolerance for the iterative method in neutrons/cm2/s = '));
        %   Checking The values
        while tol<=0 || tol>=1
        disp(sprintf('\nError!!! The tolerance is either very small, negative or very large'));
        tol=input(sprintf('Please enter the tolerance in neutrons/cm2/s again = '));
        end


        %   Entering the initial values of the flux 
        for i=1:m
            I_0(i,1)=1;
        end

        disp(sprintf('\n===========================================================================\n'));
        disp(sprintf('\t\t\t\t\t\t\t{Jacobi Iterative Method}\n'));


        %  Testing if there is iterative solution
        AbsR= abs(R);
        solution2=1;
        for i=1:m
            if AbsR(i,i) < (sum(AbsR(i,:))-AbsR(i,i))
                disp(sprintf('\n\t**Error!The diagonal magnitude in row(%d) is less the sumation of all magnitude of the row.\n\t  There''s no iterative solution.\n\t  End of the programme.',r));
                solution2=0;
                return      %   end of the programme
            end
        end

        if solution2==1
            disp(sprintf('\n**There''s an iterative solution for this circuit.\n'));
        end


        I_Jac(:,1)=I_0;

        % Initializing I(k+1) to the initial guess
        I = I_0;
        I_old = I_0;        
        % Initializing the Norm with a value greater than the tolerance (avoiding the loop to be aborted from the first iteration)
        n=100;

        % Counting the number of iterations (k)
        k = 0;


        while n >= tol

            k = k + 1;
            for i=1:m
                segma = 0;
                for j=1:m
                    if j ~= i
                        segma = segma + R(i,j)*I_old(j);
                    end
                end
                I(i) = (1/R(i,i))*(V(i)-segma);
            end

            % Calculating the infinite norm between I(k+1) & I(k)
            n = max(abs(I-I_old));    
            I_Jac(:,k+1)=I(:,1);
            I_old = I;

        end

        %  displaying the number of iterations
        disp(sprintf('\nThe number of iterations done (K): %d\n\n',k));
        disp(sprintf('The iterations results are in the table.(All the fluxes and the Norm are in (neutrons/cm2/s))\n\n'));

        for j=1:k
            Norm(j+1)=norm(I_Jac(:,j+1)-I_Jac(:,j),inf);
        end

        disp(sprintf('\nFlux and norm in the last iteration:\n'));
        Flux= I_Jac(:,k+1)
        norm= Norm(:,k+1)


        figure
        plot(0:1:k,[I_Jac ; Norm]);
        title('Jacobi Iterative Method');
        xlabel('Iterations');
        ylabel('Flux & Norm (neutrons/cm2/s)');
        grid;

    
elseif matsol==4
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%     Gauss-Seidel Iterative Method      %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %   Entering the value of the tolerance
    tol=input(sprintf('\n\nEnter the tolerance for the iterative method in neutrons/cm2/s = '));
    %   Checking The values
    while tol<=0 || tol>=1
    disp(sprintf('\nError!!! The tolerance is either very small, negative or very large'));
    tol=input(sprintf('Please enter the tolerance in neutrons/cm2/s again = '));
    end


    %   Entering the initial values of the flux 
    for i=1:m
        I_0(i,1)=1;
    end

    disp(sprintf('\n===========================================================================\n'));
    disp(sprintf('\t\t\t\t\t\t{Gauss-Seidel Iterative Method}\n'));


    
        %  Testing if there is iterative solution
    AbsR= abs(R);
    solution2=1;
    for i=1:m
        if AbsR(i,i) < (sum(AbsR(i,:))-AbsR(i,i))
            disp(sprintf('\n\t**Error!The diagonal magnitude in row(%d) is less the sumation of all magnitude of the row.\n\t  There''s no iterative solution.\n\t  End of the programme.',r));
            solution2=0;
            return      %   end of the programme
        end
    end

    if solution2==1
        disp(sprintf('\n**There''s an iterative solution for this circuit.\n'));
    end
    
    
    

    I_GS(:,1)=I_0;

    % Initializing I(k+1) to the initial guess
    I = I_0;
    I_old = I_0;        
    % Initializing the Norm with a value greater than the tolerance (avoiding the loop to be aborted from the first iteration)
    n2=100;

    % Counting the number of iterations (k)
    L = 0;


    while n2 >= tol

        L = L + 1;

        for i=1:m
            segma_1 = 0;
            for j=1:(i-1)
                segma_1 = segma_1 + R(i,j)*I(j);
            end

            segma_2 = 0;
            for j=(i+1):m
                segma_2 = segma_2 + R(i,j)*I_old(j);
            end

            I(i) = (1/R(i,i))*(V(i)-segma_1-segma_2);

        end

        % Calculating the infinite norm between I(k+1) & I(k)
        n2 = max(abs(I-I_old));    
        I_GS(:,L+1)=I(:,1);
        I_old = I;

    end

    %  Displaying the number of iterations
    disp(sprintf('\nThe number of iterations done (K): %d\n\n',L));
    disp(sprintf('The iterations results are in the table.(All the flux values and the Norm are in (neutrons/cm2/s))\n\n'));

    for j=1:L
        Norm2(j+1)=norm(I_GS(:,j+1)-I_GS(:,j),inf);
    end
    
    
    disp(sprintf('\nFlux and norm in the last iteration:\n'));
    Flux=I_GS(:,L+1)
    norm=Norm2(:,L+1)

    %  Ploting the flux & the norm

    figure
    plot(0:1:L,[I_GS ; Norm2]);
    title('Gauss-Seidel Iterative Method');
    xlabel('Iterations');
    ylabel('Flux & Norm (neutrons/cm2/s)');
    grid;

    
else matsol==5
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%     Successive Relaxation Method      %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %   Entering the value of the tolerance
    e=input(sprintf('\n\nEnter the tolerance for the iterative method in neutrons/cm2/s = '));
    %   Checking The values
    while e<=0 || e>=1
    disp(sprintf('\nError!!! The tolerance is either very small, negative or very large'));
    e=input(sprintf('Please enter the tolerance in neutrons/cm2/s again = '));
    end
    
    %  Entering the relaxation factor and sending it to the output file
    a=input(sprintf('\n\nEnter the relaxation factor (alpha) = '));
    while a<=0 || a>=2 %  checking the value
        disp(sprintf('\nError!The relaxation factor (alpha) must be between 0 and 2'));
        a=input(sprintf('Please enter the relaxation factor (alpha) again = '));
    end


    %   Entering the initial values of the flux 
    for i=1:m
        I_0(i,1)=1;
    end

    disp(sprintf('\n===========================================================================\n'));
    disp(sprintf('\t\t\t\t\t{Successive Relaxation Method}\n'));

    I_old=I_0;
    k=0;
    I_suc(:,1)=I_0;
    I_new=zeros(m,1);
    n=999;

    %  Testing if there is iterative solution
    AbsR= abs(R);
    solution2=1;
    for i=1:m
        if AbsR(i,i) < (sum(AbsR(i,:))-AbsR(i,i))
            disp(sprintf('\n\t**Error!The diagonal magnitude in row(%d) is less the sumation of all magnitude of the row.\n\t  There''s no iterative solution.\n\t  End of the programme.',r));
            solution2=0;
            return      %   end of the programme
        end
    end

    if solution2==1
        disp(sprintf('\n**There''s an iterative solution for this circuit.\n'));
    end
    
    
    

    %  The Algorithm
    while n>e

        k=k+1;

        for i=1:m

            sum1=0;
            sum2=0;

            if i>1
                for j=1:i-1
                    sum1=sum1+R(i,j)*I_new(j,1);
                end
            end

            if i<m
                for j=i+1:m
                    sum2=sum2+R(i,j)*I_old(j,1);
                end
            end

            I_new(i,1)= ((1-a)*I_old(i,1))+(a*(V(i)-sum1-sum2))/R(i,i);

        end

        n=max(abs(I_new-I_old));
        I_suc(:,k+1)=I_new(:,1);
        I_old=I_new;

    end

    %  Displaying the number of iterations
    disp(sprintf('\nThe number of iterations done (K): %d\n\n',k));
    disp(sprintf('The iterations results in the table.(All the currents and the Norm are in (neutrons/cm2/sec))\n\n'));

    for j=1:k
        Norm(j+1)=norm(I_suc(:,j+1)-I_suc(:,j),inf);
    end

    disp(sprintf('\nFlux and norm in the last iteration:\n'));
    Flux=I_suc(:,k+1)
    norm=Norm(:,k+1)

    %  Ploting the flux & the norm

    plot(0:1:k,[I_suc ; Norm]);
    title('Successive Relaxation Plot');
    xlabel('Iterations');
    ylabel('Flux & Norm (neutrons/cm2/s)');
    grid;

end

conex=input(sprintf('\n\nPress 1 and enter to continue. \nPress 0 and Enter to stop the program. \n1 or 0? '));

end

%  Massage for the user
disp(sprintf('\n\n\t\t\t\t\t***END OF THE PROGRAMME***'));

