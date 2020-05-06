clc,clear,close all
warning off
format longG
% 万有引力定律算法
% GSA 参数
G0 = 100;
a = 20;
Vmin = -1;
Vmax = 1;
maxiter = 200;  % 迭代次数
sizepop = 20;  % 种群数量
popmin1 = -1;  popmax1 = 1; % x1
popmin2 = -1;  popmax2 = 1; % x2
% 初始化种群
for i=1:sizepop
    x1 = popmin1 + (popmax1-popmin1)*rand;
    x2 = popmin2 + (popmax2-popmin2)*rand;
    pop(i,1) = x1;
    pop(i,2) = x2;
    fitness(i) = fun([x1,x2]);
end
% 记录一组最优值
[bestfitness,bestindex]=min(fitness);
zbest=pop(bestindex,:);   % 全局最佳
gbest=pop;                % 个体最佳
fitnessgbest=fitness;     % 个体最佳适应度值
fitnesszbest=bestfitness; % 全局最佳适应度值
V = [0,0];
% 迭代寻优
for i=1:maxiter
   
    Gt = G0*exp(-a*i/maxiter);        % 万有引力定律常数
   
    [mfitness,index] = sort(fitness); % 从小到大排序
    for j=1:sizepop
        mt(j) = (fitness(j)-mfitness(end))./(mfitness(1)-mfitness(end)+eps);
    end
   
    for j=1:sizepop
        Mt = abs(mt(j))./sum(abs(mt))+eps;
        Rij = norm( pop(j,:)-gbest(j,:) ) +eps;
        Fij = rand(1,2).*( Gt*Mt*Mt/(Rij+eps) .*(pop(j,:)-gbest(j,:)+eps ) );
        Va = Fij./Mt;
        
        % 速度更新
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
        
        % 位置更新
        pop(j,:) = pop(j,:) + V;
        % x1  越界限制
        if pop(j,1)>popmax1
            pop(j,1)=popmax1;
        end
        if pop(j,1)<popmin1
            pop(j,1)=popmin1;
        end
        % x2  越界限制
        if pop(j,2)>popmax2
            pop(j,2)=popmax2;
        end
        if pop(j,2)<popmin2
            pop(j,2)=popmin2;
        end
        
        % 适应度更新
        fitness(j) = fun(pop(j,:));
        
        % 比较  个体间比较
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
disp('最优解')
disp(zbest)
fprintf('\n')

figure('color',[1,1,1])
plot(fitness_iter,'ro-','linewidth',2)

figure('color',[1,1,1])
loglog(fitness_iter,'ro-','linewidth',2)
axis tight
grid on