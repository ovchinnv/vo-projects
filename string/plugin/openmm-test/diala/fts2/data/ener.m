
styles={'r','g','b','m','c','k','y'};

ind=[10:10:50];

figure; hold on; box on; grid on;

leg={};
for i=1:length(ind)
 e=load(['string',num2str(ind(i)),'.dat']);
 rep=e(:,1);
 tote=e(:,2);
 plot(rep,tote, char(styles(mod(i-1,length(styles))+1)),'linewidth',1.5);
 leg=[leg, {['iteration ', num2str(ind(i))] }];
end 
 
legend(leg,-1);
 
xlabel('\it Replica', 'fontsize',14);
ylabel('\it Potential energy (kcal/mol)','fontsize',14);

set(gcf,'Paperpositionmode','auto');
print(gcf, '-depsc2','zts.eps');
%print(gcf, '-djpeg100','zts.jpg');
