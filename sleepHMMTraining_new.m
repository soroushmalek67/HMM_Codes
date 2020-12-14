function [HMM,EMIS,TRANS,Likelihood]=sleepHMMTraining_new(Obs,States,numNeurons,Repetitions,Iterations)

   
   EMIS       = cell(1,Repetitions);
   TRANS      = cell(1,Repetitions);
   Likelihood = zeros(1,Repetitions);
   
      for i=1:Repetitions
        
            display(['repetition : ' num2str(i)])
            
        % Transition matrix initialization:
        
        D=0.99+0.009*rand(States,1);
        TRANS_GUESS = diag(D);
        for s=1:States
        TRANS_GUESS(s,setxor(1:States,s)) = (1-D(s))/(States-1);    
        end        
    
         % Emission matrix initialization:
         
        EMIS_GUESS = rand(States,numNeurons+1);
        EMIS_GUESS(1,:) = 0.002*rand(1,numNeurons+1);
        EMIS_GUESS(1,1) = .99; 

        [TRANS_EST, EMIS_EST] = hmmtrain(Obs, TRANS_GUESS, EMIS_GUESS,...
            'Algorithm', 'BaumWelch','MAXITERATIONS',Iterations);
        
          if isnumeric(Obs)
          % Calculate likelihoods
            [P, L] = hmmdecode(Obs, TRANS_EST, EMIS_EST);
            Likelihood(i)=L;
            EMIS{i}=EMIS_EST;
            TRANS{i}=TRANS_EST;
          else
            for k=1:length(Obs)
            [gP, gL] = hmmdecode(Obs{k}, TRANS_EST, EMIS_EST);
            L(k)=gL;
            end
            Likelihood(i)=sum(L);
            EMIS{i}=EMIS_EST;
            TRANS{i}=TRANS_EST;            
          end
          
        
      end
      
%         ML=min(Likelihood);
        ML=max(Likelihood);
            
        HMM.EMISSION=EMIS{Likelihood==ML};
        HMM.TRANSITION=TRANS{Likelihood==ML};
    
        
return