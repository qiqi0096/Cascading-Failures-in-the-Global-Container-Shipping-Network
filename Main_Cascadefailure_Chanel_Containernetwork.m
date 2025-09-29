clc
clear all
tic

load finaltable1_L_noweight
load finaltable1_L_weight
load finaltable1_P_noweight
load portname

%% simulation model
x_firstload=1; 
x_attract=2; 
x_capacity1_min=0.1; 
x_capacity1_max=1;
x_capacity2_min=0;
x_capacity2_max=1.5; 
stepx=0.1;
a=(x_capacity2_max-x_capacity2_min)/stepx+1;
b=(x_capacity1_max-x_capacity1_min)/stepx+1;

simulationdata1=zeros(a,b); 
simulationdata2=zeros(a,b); 
step1=zeros(a*b,1);
step2=zeros(a*b,1);

x_capacity1_sum=[];
x_capacity2_sum=[];
A=finaltable1_P_noweight(2:end,2:end); % Space P
A_weightroute=finaltable1_L_weight(2:end,2:end); % weight
A_spaceL=finaltable1_L_noweight(2:end,2:end); % Space L
Num_Nodes=size(A,1);
port=(1:Num_Nodes)';
[degree_all,Ave_degree00,Node_RATE00,Paths,Shortest_L00,Eff_I00,Diameter00,Giant_tu00]=Calculatefunction_Direct(A,port,port);% network metrics
x_capacity2=x_capacity2_min;
for m=1:a
    % m=1;n=1;
    x_capacity1=x_capacity1_min; 
    for n=1:b
        F=x_firstload*degree_all; % F
        C=F+x_capacity1*(F.^x_capacity2);% C
        test=F./C;
        numbertext1=sum(sum(test<=0.68));
        numbertext2=sum(sum(test>=0.72));
        if numbertext1==0 && numbertext2==0
            x_capacity1_sum=[x_capacity1_sum;x_capacity1];
            x_capacity2_sum=[x_capacity2_sum;x_capacity2];
        end
        x_capacity1=x_capacity1+stepx; 
    end
    x_capacity2=x_capacity2+stepx; 
end
x_capacity2=x_capacity2_min; 
for m=1:a
    % m=1;n=1;
    x_capacity1=x_capacity1_min;
    for n=1:b
        A=finaltable1_P_noweight(2:end,2:end); % Space P
        A_weightroute=finaltable1_L_weight(2:end,2:end); % weight
        A_spaceL=finaltable1_L_noweight(2:end,2:end); % Space L
        Num_Nodes=size(A,1);
        portreal=(1:Num_Nodes)'; 
        port=(1:Num_Nodes)';

        %% Risk1 
        risk=1;
        Delete_allports=[276 277 278 279 280 283 284 285 306 307 308 309 310 312 820 823 824 825 841 864 865]; % Delete ports
        Delete_allports(Delete_allports>=432)=Delete_allports(Delete_allports>=432)-1; 
        De_choose=99;
        alt_port=[835 836 837 838 842 843 844 845 846 847 848 849 850 851 852 853 854 855 856 857 858 859 860 861 862 863 866 867 868 869 870 871 872 873 874 877 878 879 880 881 882 883 884 885 886 887 888 889 890 891 892 893 900 901 902 903 904 905 906];
        alt_port(alt_port>=432)=alt_port(alt_port>=432)-1; 
        % Alternative Ports: Initial sequence numbers for ports along the southern African coastline stretching from Kenya in the east to Morocco in the west, 
        % as well as ports in island nations such as Madagascar. (If the results are unreasonable, consider adjusting parameters, as Africa's infrastructure is
        % less developed than other regions and port resources cannot be effectively mobilized.)

        %%  Risk2 
        % risk=2;
        % Delete_allports=[276 277 278 279 280 283 284 285 306 307 308 309 310 312 820 823 824 825 841 864 865 286 287 288 289 290 291 292 293 294 295 296 297 298 299 300 301 302 303 304 305 313 314 315 316 816 817 818 819 821 822 875 876]; % 失效红海+东地中海海域（利比亚、希腊以东，目的是区分欧洲内部港口之间仍然可以贸易往来，但是由西向东穿越地中海不行）港口的索引号（直接删除）
        % Delete_allports(Delete_allports>=432)=Delete_allports(Delete_allports>=432)-1; 
        % De_choose=99; 
        % alt_port=[835 836 837 838 842 843 844 845 846 847 848 849 850 851 852 853 854 855 856 857 858 859 860 861 862 863 866 867 868 869 870 871 872 873 874 877 878 879 880 881 882 883 884 885 886 887 888 889 890 891 892 893 900 901 902 903 904 905 906]; % 替代港口的初始序号
        % alt_port(alt_port>=432)=alt_port(alt_port>=432)-1; 
        % % Alternative Ports: Initial sequence numbers for ports along the southern African coastline stretching from Kenya in the east to Morocco in the west, 
        % % as well as ports in island nations such as Madagascar. (If the results are unreasonable, consider adjusting parameters, as Africa's infrastructure is
        % % less developed than other regions and port resources cannot be effectively mobilized.)
        %%  Risk3 
        % risk=3;
        % Delete_allports=[276 277 278 279 280 283 284 285 306 307 308 309 310 311 312 820 823 824 825 841 864 865 286 287 288 289 290 291 292 293 294 295 296 297 298 299 300 301 302 303 304 305 313 314 315 316 816 817 818 819 821 822 875 876]; % 失效红海+东地中海海域（利比亚、希腊以东，目的是区分欧洲内部港口之间仍然可以贸易往来，但是由西向东穿越地中海不行）港口的索引号（直接删除）
        % Delete_allports(Delete_allports>=432)=Delete_allports(Delete_allports>=432)-1; 
        % De_choose=99; 
        % alt_port=[835 836 837 838 842 843 844 845 846 847 848 849 850 851 852 853 854 855 856 857 858 859 860 861 862 863 866 867 868 869 870 871 872 873 874 877 878 879 880 881 882 883 884 885 886 887 888 889 890 891 892 893 900 901 902 903 904 905 906]; % 替代港口的初始序号
        % alt_port(alt_port>=432)=alt_port(alt_port>=432)-1;
        [degree_all,Ave_degree00,Node_RATE00,Paths,Shortest_L00,Eff_I00,Diameter00,Giant_tu00]=Calculatefunction_Direct(A,port,port);

        %% F
        F=x_firstload*degree_all;

        %% C
        C=F+x_capacity1*(F.^x_capacity2);

        %% Strength
        Path1=Paths;
        idy= Path1==0;
        Path1(idy)=Inf;
        Strength=x_attract*A_weightroute'.*(C-F)./Path1';
        Strength=Strength';

        %% f
        Strength=Strength.*A_spaceL; 
        f1=Strength./sum(Strength,2);
        idy= isnan(f1);
        f1(idy)=0;
        if De_choose==99 
            f2=zeros(1,Num_Nodes);
            f2(alt_port)=(degree_all(alt_port,:)./sum(degree_all(alt_port),1))'; 
            f2=ones(Num_Nodes,Num_Nodes).*f2;
            idy= isnan(f2);
            f2(idy)=0;
        end

        %% Initial attack node removal
        Fail_All=[]; % Fail_Allport
        Remove_Free=[]; 
        Nodes_allRemove=[]; 
        Main_Nodes=1:Num_Nodes; 
        NumberOfNodes_Temp=Num_Nodes; 
        
        %% attack
        F00=F;
        C00=C;
        f00=f1;
        A00=A;
        Fail_Sum=0;
        port00=port;
        Nodes_allRemove=[Nodes_allRemove;Delete_allports']; 

        port00(Delete_allports,:)=[];  
        A00=A(port00,port00);  
        [~,~,Alone_nodes_index,Alone_nodes_Num]=Degree_Distribution_Nofigure(A00);  

        if ~isempty(Alone_nodes_index)==1 % 
            Fail_Sum=Fail_Sum+Alone_nodes_Num; % 
            Alone_nodes=port00(Alone_nodes_index,:);  
            Nodes_allRemove=[Nodes_allRemove;Alone_nodes]; 
            port00(Alone_nodes_index,:)=[]; 
            A00=A(port00,port00);  
        end

        if De_choose==1
            %% Calculate load redistribution caused by failed nodes (removed nodes + isolated nodes)
            f00=f1(Nodes_allRemove,:); 
            f00(:,Nodes_allRemove)=[]; 
            f00=f00./sum(f00,2);
            idz= isnan(f00);
            f00(idz)=0;
            F_Remove=F00(Nodes_allRemove,:);  
            F00(Nodes_allRemove,:)=[];  
            F00=F00+(F_Remove'*f00)';  
            C00=C(port00,:);
     
            NodesFailure=F00>C00; 

        else % Select alternative ports of call for substitute routes

            f00=f2(Nodes_allRemove,:);  
            f00(:,Nodes_allRemove)=[]; 
            f00=f00./sum(f00,2); 
            idz= isnan(f00);
            f00(idz)=0;
            F_Remove=F00(Nodes_allRemove,:);  
            F00(Nodes_allRemove,:)=[];  
            F00=F00+(F_Remove'*f00)'; 
            C00=C(port00,:); 

            NodesFailure=F00>C00;  

        end

        %% In the event of a cascading failure, remove the affected node and transfer the load exceeding its capacity; if an isolated node occurs, transfer all its load.
        % Whether bypassing ports or using alternative ports, load distribution during cascading failures follows the f1 gravity model.
        while sum(NodesFailure)>0  
            step1(b*(m-1)+n)=step1(b*(m-1)+n)+1; 
            Casfail_Index=find(NodesFailure); 
            Fail_Sum=Fail_Sum+numel(Casfail_Index); 
            Nodes_Casfail = port00(Casfail_Index,:); 
            Nodes_allRemove=[Nodes_allRemove;Nodes_Casfail];
            New_Remove=Nodes_Casfail;
            port00(Casfail_Index,:)=[]; 
            % F_Remove=F00(Casfail_Index,:);  
            F_Remove=(F00(Casfail_Index,:)-C00(Casfail_Index,:)); 
            F00(Casfail_Index,:)=[]; 
            A00=A(port00,port00); 
            [~,~,Alone_nodes_index,Alone_nodes_Num]=Degree_Distribution_Nofigure(A00);  

            if ~isempty(Alone_nodes_index)==1 
                Fail_Sum=Fail_Sum+Alone_nodes_Num; 
                Alone_nodes=port00(Alone_nodes_index,:);  
                Nodes_allRemove=[Nodes_allRemove;Alone_nodes]; 
                New_Remove=[Nodes_Casfail;Alone_nodes];
                port00(Alone_nodes_index,:)=[]; 
                F_Remove=[F_Remove;F00(Alone_nodes_index,:)]; 
                F00(Alone_nodes_index,:)=[];  
                A00=A(port00,port00); 
            end

            %% Calculate load redistribution caused by failed nodes (cascade failures + isolated nodes): 
            % Note that load F must be updated based on the previous iteration.
            f00=f1(New_Remove,:); 
            f00(:,Nodes_allRemove)=[]; 
            f00=f00./sum(f00,2); 
            idz= isnan(f00);
            f00(idz)=0;
            F00=F00+(F_Remove'*f00)'; 
            C00=C(port00,:); 

            NodesFailure=F00>C00; 
        end

        %% 
        if risk==3
            De_choose=1; 
            Delete_allports=[317 318 319 320 321 322 323 325 326 330 331 332 333 334 335 336]; 
            Delete_allports(Delete_allports>=432)=Delete_allports(Delete_allports>=432)-1; 
            [Sameports1,index1,~]=intersect(Delete_allports',Nodes_allRemove);
            Delete_allports(index1)=[];
            [Sameports2,~,index2]=intersect(Delete_allports',port00);
            
            Nodes_allRemove=[Nodes_allRemove;Delete_allports'];  
            New_Remove=Delete_allports';
            port00(index2,:)=[];  

            F_Remove=F00(index2,:); 
            F00(index2,:)=[]; 
            A00=A(port00,port00); 
            [~,~,Alone_nodes_index,Alone_nodes_Num]=Degree_Distribution_Nofigure(A00); 

            if ~isempty(Alone_nodes_index)==1 
                Fail_Sum=Fail_Sum+Alone_nodes_Num;
                Alone_nodes=port00(Alone_nodes_index,:); 
                Nodes_allRemove=[Nodes_allRemove;Alone_nodes]; 
                New_Remove=[Delete_allports';Alone_nodes];
                port00(Alone_nodes_index,:)=[]; 
                F_Remove=[F_Remove;F00(Alone_nodes_index,:)]; 
                F00(Alone_nodes_index,:)=[];  
                A00=A(port00,port00); 
            end

            %% Calculate load redistribution caused by failed nodes (cascade failures + isolated nodes): 
            % Note that load F must be updated based on the previous iteration.
            f00=f1(New_Remove,:);   
            f00(:,Nodes_allRemove)=[]; 
            f00=f00./sum(f00,2);
            idz= isnan(f00);
            f00(idz)=0;
            F00=F00+(F_Remove'*f00)'; 
            C00=C(port00,:);  

           
            NodesFailure=F00>C00;  

            %%   In the event of a cascading failure, remove the affected node and transfer the load exceeding its capacity; 
            % should an isolated node occur, transfer all its load.

            while sum(NodesFailure)>0  
                step2(b*(m-1)+n)=step2(b*(m-1)+n)+1; 
                Casfail_Index=find(NodesFailure);
                Fail_Sum=Fail_Sum+numel(Casfail_Index); 
                Nodes_Casfail = port00(Casfail_Index,:); 
                Nodes_allRemove=[Nodes_allRemove;Nodes_Casfail]; 
                New_Remove=Nodes_Casfail;
                port00(Casfail_Index,:)=[];

                F_Remove=(F00(Casfail_Index,:)-C00(Casfail_Index,:)); 
                F00(Casfail_Index,:)=[];  
                A00=A(port00,port00); 
                [~,~,Alone_nodes_index,Alone_nodes_Num]=Degree_Distribution_Nofigure(A00);

                %%   Locate isolated nodes resulting from the removal of cascading failure nodes
                if ~isempty(Alone_nodes_index)==1 
                    Fail_Sum=Fail_Sum+Alone_nodes_Num;
                    Alone_nodes=port00(Alone_nodes_index,:); 
                    Nodes_allRemove=[Nodes_allRemove;Alone_nodes];
                    New_Remove=[Nodes_Casfail;Alone_nodes];
                    port00(Alone_nodes_index,:)=[];
                    F_Remove=[F_Remove;F00(Alone_nodes_index,:)]; 
                    F00(Alone_nodes_index,:)=[]; 
                    A00=A(port00,port00);
                end

                %% In the event of a cascading failure, remove the affected node and transfer the load exceeding its capacity; if an isolated node occurs, transfer all its load.
                % Whether bypassing ports or using alternative ports, load distribution during cascading failures follows the f1 gravity model.
                f00=f1(New_Remove,:);  
                f00(:,Nodes_allRemove)=[];  
                f00=f00./sum(f00,2); 
                idz= isnan(f00);
                f00(idz)=0;
                F00=F00+(F_Remove'*f00)';
                C00=C(port00,:); 

                NodesFailure=F00>C00; 
            end
        end

        [degree_all_after,Ave_degree00_after,Node_RATE00_after,Paths_after,Shortest_L00_after,Eff_I00_after,Diameter00_after,Giant_tu00_after]=Calculatefunction_Direct(A00,port00,port);
        simulationdata1(m,n)=Giant_tu00_after;
        simulationdata2(m,n)=Eff_I00_after;
        x_capacity1=x_capacity1+stepx; 
    end
    x_capacity2=x_capacity2+stepx; 
end
% save Risk2
toc