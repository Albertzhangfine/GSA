clc,clear,close all
warning off
format longG
% �������������㷨
% GSA ����
G0 = 100;
a = 20;
Vmin = -1;
Vmax = 1;
maxiter = 200;  % ��������
sizepop = 20;  % ��Ⱥ����
popmin1 = -1;  popmax1 = 1; % x1
popmin2 = -1;  popmax2 = 1; % x2
% ��ʼ����Ⱥ
for i=1:sizepop
    x1 = popmin1 + (popmax1-popmin1)*rand;
    x2 = popmin2 + (popmax2-popmin2)*rand;
    pop(i,1) = x1;
    pop(i,2) = x2;
    fitness(i) = fun([x1,x2]);
end
% ��¼һ������ֵ
[bestfitness,bestindex]=min(fitness);
zbest=pop(bestindex,:);   % ȫ�����
gbest=pop;                % �������
fitnessgbest=fitness;     % ���������Ӧ��ֵ
fitnesszbest=bestfitness; % ȫ�������Ӧ��ֵ
V = [0,0];
% ����Ѱ��
for i=1:maxiter
   
    Gt = G0*exp(-a*i/maxiter);        % �����������ɳ���
   
    [mfitness,index] = sort(fitness); % ��С��������
    for j=1:sizepop
        mt(j) = (fitness(j)-mfitness(end))./(mfitness(1)-mfitness(end)+eps);
    end
   
    for j=1:sizepop
        Mt = abs(mt(j))./sum(abs(mt))+eps;
        Rij = norm( pop(j,:)-gbest(j,:) ) +eps;
        Fij = rand(1,2).*( Gt*Mt*Mt/(Rij+eps) .*(pop(j,:)-gbest(j,:)+eps ) );
        Va = Fij./Mt;
        
        % �ٶȸ���
        V = rand(1,2).*V - Va;
        % V--x1
        if V(1,1)>Vmax
            V(1,1)=Vmax;
        end
        if V(1,1)<Vmin
            V(1,1)=Vmin;
        end
        % V--x2
        if V(1,2)>Vmax
            V(1,2)=Vmax;
        end
        if V(1,2)<Vmin
            V(1,2)=Vmin;
        end
        
        % λ�ø���
        pop(j,:) = pop(j,:) + V;
        % x1  Խ������
        if pop(j,1)>popmax1
            pop(j,1)=popmax1;
        end
        if pop(j,1)<popmin1
            pop(j,1)=popmin1;
        end
        % x2  Խ������
        if pop(j,2)>popmax2
            pop(j,2)=popmax2;
        end
        if pop(j,2)<popmin2
            pop(j,2)=popmin2;
        end
        
        % ��Ӧ�ȸ���
        fitness(j) = fun(pop(j,:));
        
        % �Ƚ�  �����Ƚ�
        if fitness(j)<fitnessgbest(j)
            fitnessgbest(j) = fitness(j);
            gbest(j,:) = pop(j,:);
        end
        if fitness(j)<bestfitness
            bestfitness = fitness(j);
            zbest =  pop(j,:);
        end
    end
    fitness_iter(i) = bestfitness;
end
disp('���Ž�')
disp(zbest)
fprintf('\n')

figure('color',[1,1,1])
plot(fitness_iter,'ro-','linewidth',2)

figure('color',[1,1,1])
loglog(fitness_iter,'ro-','linewidth',2)
axis tight
grid on