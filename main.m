%% The script corresponds to the Algorithm 5 with r = ( 3 1) 

clc;
clear ;
close all ;

%% The Results.txt
% This file contains the ultimate SDR values for different channel parameters delta and epsilon. 
% Also, since the proposed ACOSQ is an iterative algorithm, the distirtion
% value for every iteration is provided in this file for given a epsilon and delta. 

FileID = fopen ('Results.txt' , 'a') ;
alpha = 500000 ;



%% Number of Quantization level
numLevel = 16 ;

%% Channel's cross-over probability epsilon
epsilon = unique ([10^-5 : 2 * 10^-5 : 10^ -4 , 10^-4 : 10^ -4 : 10^ -3 , 10^-3 ,0.005 0.01  0.05  0.1]);

% Since the design converges to a locally optimal solution, to avoid
% bad local optimums, we use a increase-decrease method. 
SIZE = length(epsilon) ;
noise = [1 : SIZE , SIZE : -1 : 1 , 1 : SIZE ] ;


%% Distortion parameter
D_4 = zeros(length(noise) , 1) ;

%% SDR parameters
Final_SDR_2 = zeros(24 , length(noise)) ;
SDR_2 = zeros(length(noise) , 1) ;
SDR_1 = zeros(length(noise) , 1) ;

%% Noise correlation
% The variable delta determines the amount of noise correlation. 
for delta = [0 5 10]
    % Set the channel's cross-over probability
    for k = SIZE : length(noise)
        i = noise(k);
        
        Pr_1 = [1 - epsilon(i) , epsilon(i) ;
            epsilon(i) , 1 - epsilon(i)] ;
        
        Pr_z = [(1 - epsilon(i) + delta) / (1 + delta)  , epsilon(i) / (1 + delta) ;
            (1 - epsilon(i)) / (1 + delta)  , (epsilon(i) + delta) / (1 + delta)] ;
        
        % Find the channels transition distribution for a given number of
        % quantization levels numLevel
        Pr = Channel_with_Memory(numLevel , epsilon(i) , delta) ;
        
        % Set up the parameters for computing the integrals. Note that in
        % this script all integrals are computed numerically using the
        % Riemann summation. 
        delta_u = 8 / 2^ 11 ;
        T_1(: , 1) = -4 : delta_u : 4 ;
        u = T_1(: , 1) ;
        
        % Compute the source pdf. We herein consider a zero-mean
        % unit-variance Gaussian source distribution. 
        
        f =  1 ./ (sqrt (2 .* pi)) .* exp (-u .^ 2 ./ 2) ;
        f = f ./ (sum(f) .* delta_u) ;
        counter = 0 ;
        
        % As noted in the Thesis, the ultimate codebook otabined in the
        % last step of the ACOSQ decribed in Section 4.1 is used as the 
        % initial state of the proposed ACOSQ with identical noise
        % correlation and smallest cross-over probability. 
        if (k == 1 )
            LOAD = ['ACOSQ_3_1_codebook_delta_' num2str(delta)] ;
            load (LOAD) ;
            codebook_1 = hold_codebook_4 ;
        else
            % We slightly increase the channel's cross-over probability,
            % setting the codebook from the system with small epsilon as
            % the initial state of the system with new epsilon. 
            load codebook_1
        end
        
        % The first step the proposed ACOSQ. Design a 4 bit COSQ 
        [SDR_1(k) , ~ , T_1 , codebook_1] = ACOSQ_step_1(f , Pr , numLevel , T_1 , codebook_1 , delta_u , 1 : 16) ;
        
        % save the codebook to initialize the system with the next value of
        % epsilon. 
        save('codebook_1'  , 'codebook_1') ;
        
        % save the partition set and codebook for the computing the experimental results.  
        Data = ['T\T_1_k_' num2str(k) '_delta_' num2str(delta)] ;
        save(Data , 'T_1' , 'codebook_1') ;
        
        % Exhaustively search for the 3 bits to transmit over the channel. 
        for bit_pattern = [1 2 3 ; 1 2 4 ; 1 3 4 ; 2 3 4]'
            permutations = perms(bit_pattern) ;
            for permute_index= 1 : length(permutations)
                permuted_sequence = permutations(permute_index , :) ;
                counter = counter + 1 ;
                
                % Scramble the partition indexes based on the three bits chosen for
                % tranmision such that the bits chosen for transmission are
                % always located first three places in the channel input sequence.  
                [codebook_2 , init_T_1] = scramble_labels_step_1(f , Pr , T_1 , numLevel , permuted_sequence) ;
      
                
                % The second step of the proposed ACOSQ design. In this step,
                % the remaining bits i.e., the fourth bit is generated adaptive to
                % the received bits y_1y_2y_3 corresponding to the previous
                % transmision in step one.
            
                [SDR_2(k) , T_2 , codebook_2] = ACOSQ_step_2(f , Pr , Pr_z , codebook_2 , init_T_1 , numLevel , delta_u) ;
                Final_SDR_2(counter , k) = SDR_2(k) ;
                
                % save the partition set and codebook for the computing the experimental results.  
                Data = ['T\T_2_k_' num2str(k) '_delta_' num2str(delta) , '_counter_' num2str(counter)] ;
                save (Data , 'T_2' , 'codebook_2') ;
            end
            fprintf (FileID , 'i = %d\n' , i) ;
        end
    end
    myFinal_SDR_2 = max(Final_SDR_2 , [] , 1) ;
    % Pick the best SDR value after the end of the so-called
    % increase-decrease method. 
    final_SDR_2 = zeros(SIZE , 1) ;
    for i = 1 : SIZE
        index = find (noise == i) ;
        hold_var  = myFinal_SDR_2(index) ;
        hold_var = hold_var (:) ;
        [final_SDR_2(i)  , index] = max(hold_var) ;
        fprintf (FileID , 'i = %d\n' , i) ;
        fprintf (FileID , '\nfinal_SDR_2 = %f' , final_SDR_2(i)) ;
        fprintf (FileID , '\nindex = %d\n' , index) ;
    end
    % Save the best SDR values for every channel parameters. 
    Data = ['Proposed_ACOSQ_3_1_delta_' num2str(delta)] ;
    save (Data , 'final_SDR_2' , 'myFinal_SDR_2' , 'Final_SDR_2' , 'epsilon') ;
    
end