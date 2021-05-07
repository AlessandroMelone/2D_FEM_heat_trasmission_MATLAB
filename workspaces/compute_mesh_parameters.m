%%
addpath('workspaces');
load domain_toolbox_4;
%%
total_sum = 0; %total sum of segment lenght
side_indx=[1 2;2 3;3 1];
already_consid = zeros(1,2);
number_segment = 0;
for i_t = 1:size(t,2) %for each triangle sum the lenght three segments
    for n_side = 1:3 %for each segment
        i = t(side_indx(n_side,1),i_t); %index first node of the segment
        j = t(side_indx(n_side,2),i_t); %index second node of the segment
       
        A1 = [i,j] == already_consid;
        A2 = [j,i] == already_consid;
        if(~any(or(A1,A2)))
            squares = (p(:,i)-p(:,j)).^2;
            total_sum = total_sum + sqrt(squares(1)+squares(2));
            already_consid(number_segment+1,:) = [i,j];
            number_segment = number_segment+1;
        end
    end
    
end
average = total_sum/number_segment;

% %% Result:
% disp(['Domain1:  num triangle=',num2str(size(t,2)),' num nodes=',num2str(size(p,2))]);
% disp(['          num segment=',num2str(number_segment),' average lenght=',num2str(average)]);

% disp(['Domain2:  num triangle=',num2str(size(t,2)),' num nodes=',num2str(size(p,2))]);
% disp(['          num segment=',num2str(number_segment),' average lenght=',num2str(average)]);

% disp(['Domain3:  num triangle=',num2str(size(t,2)),' num nodes=',num2str(size(p,2))]);
% disp(['          num segment=',num2str(number_segment),' average lenght=',num2str(average)]);

disp(['Domain4:  num triangle=',num2str(size(t,2)),' num nodes=',num2str(size(p,2))]);
disp(['          num segment=',num2str(number_segment),' average lenght=',num2str(average)]);

% Domain1:  num triangle=262 num nodes=150
%           num segment=71 average lenght=0.3288
% 
% Domain2:  num triangle=1048 num nodes=561
%           num segment=249 average lenght=0.16528
% 
% Domain3:  num triangle=4192 num nodes=2169
%           num segment=981 average lenght=0.086124
% 
% Domain4:  num triangle=16768 num nodes=8529
%           num segment=3734 average lenght=0.042625

