% #108 from the microglia genes is deleted.
microglia_gene=all_raw(logical(celltype(:,2)),:);
norm_microglia_gene=microglia_gene./(max(microglia_gene,[],2));

delete = [];
for i = 1:size(norm_microglia_gene,1)
    temp = min(norm_microglia_gene(i,:));
    if temp == 1
        delete=[delete,i];
    end
end
norm_microglia_gene(delete,:)=[];

microglia_adj=[];
for i = 1:size(norm_microglia_gene,1)
    for j = 1:size(norm_microglia_gene,1)
        temp_cor=corrcoef(norm_microglia_gene(i,:),norm_microglia_gene(j,:));
        microglia_adj(i,j)=temp_cor(1,2);
    end
end

%define community based on Louvain method
A=abs(microglia_adj);
gamma = 1;
k = full(sum(A));
twom = sum(k);
B = full(A - gamma*k'*k/twom);
[S,Q] = genlouvain(B);
com_color=zeros(234,3);
for i = 1:234
    com_color(i,S(i))=1;
end

%plot graph
for i = 1:size(norm_microglia_gene,1)
    microglia_adj(i,i)=0;
end
microglia_th=abs(microglia_adj)>0.31;
G = graph(microglia_th);
p=plot(G,'Layout','force');
p.NodeColor=com_color;

%plot schemaball
[S_sorted,S_idx]=sort(S);
norm_microglia_rearrange=norm_microglia_gene(S_idx,:);
microglia_name=names_sort(logical(celltype(:,2)));
microglia_name(delete)=[];
name_microglia_rearrange=microglia_name(S_idx);
microglia_adj_rearrange=[];
for i = 1:size(norm_microglia_rearrange,1)
    for j = 1:size(norm_microglia_rearrange,1)
        temp_cor=corrcoef(norm_microglia_rearrange(i,:),norm_microglia_rearrange(j,:));
        microglia_adj_rearrange(i,j)=temp_cor(1,2);
    end
end
microglia_adj_rearrange = abs(microglia_adj_rearrange);
nodeColor=zeros(length(S),3);
for i = 1:length(S)
    nodeColor(i,S_sorted(i))=1;
end
[~,rank]=sort(sum(microglia_adj_rearrange));
H=schemaball(microglia_adj_rearrange,name_microglia_rearrange,[1,1,1],[]);
H.s.CData=nodeColor(transpose(rank),:);
H.s.MarkerEdgeAlpha=1;
H.s.MarkerEdgeColor=[0.1,0.1,0.1];
H.s.LineWidth=0.5;
H.s.MarkerFaceAlpha=0.5;
idx=1:40;
idx_p=idx(~isnan(H.l));

for i = 1:15
    set(H.l(idx_p(i)),'Color',[0.1,0,0]);
    set(H.l(idx_p(i)),'LineWidth',0.1);
end
for i = 16:20
    set(H.l(idx_p(i)),'Color',[0.04*i,0.3,0.3]);
    set(H.l(idx_p(i)),'LineWidth',0.1+0.3*(i-16));
end

%Original size
H.s.SizeData=500.^0.5;

%Size by TG vs WT
TG_vs_WT=mean(microglia_gene(:,21:30),2)./mean(microglia_gene(:,53:58),2);
TG_vs_WT(delete)=[];
TG_vs_WT_rearrange=TG_vs_WT(S_idx);
TG_vs_WT_rearrange_scale=1./(1 + exp(-5.*(TG_vs_WT_rearrange-1)));
H.s.SizeData=500.^(TG_vs_WT_rearrange_scale(rank));

%Size by PLS vs TG
PLX_vs_TG=mean(microglia_gene(:,41:46),2)./mean(microglia_gene(:,53:58),2);
PLX_vs_TG(delete)=[];
PLX_vs_TG_rearrange=PLX_vs_TG(S_idx);
PLX_vs_TG_rearrange_scale=1./(1 + exp(-5.*(PLX_vs_TG_rearrange-1)));
H.s.SizeData=500.^(PLX_vs_TG_rearrange_scale(rank));


%subnetwork
sub_group2 = microglia_adj_rearrange(92:185,92:185);
for i = 1:size(sub_group2,1)
    sub_group2(i,i)=0;
end
persist_name=name_microglia_rearrange(92:185);
persist_th=abs(sub_group2)>0.7;
PERSIS=graph(persist_th);
per=plot(PERSIS,'Layout','force','NodeLabel',persist_name);

gamma = 1;
k = full(sum(persist_th));
twom = sum(k);
B_P = full(persist_th - gamma*k'*k/twom);
[S_P,Q] = genlouvain(B_P);
[~,S_p_idx]=sort(S_P);

pmap=colormap(jet);
pmap=imresize(pmap,[32,3]);
persist_color2=zeros(94,3);
for i = 1:94
    if S_P(i)==2
        persist_color2(i,:)=[1,0,0];
    elseif S_P(i)==4
        persist_color2(i,:)=[0,1,0];
    elseif S_P(i)==7
        persist_color2(i,:)=[0,0,1];
    end
end

%% excluding the PLX groups
groups_noPLX=transpose(([21:40,53:64]));
norm_microglia_gene_noPLX=norm_microglia_gene(:,groups_noPLX);
microglia_adj_noPLX=[];
for i = 1:size(norm_microglia_gene_noPLX,1)
    for j = 1:size(norm_microglia_gene_noPLX,1)
        temp_cor=corrcoef(norm_microglia_gene_noPLX(i,:),norm_microglia_gene_noPLX(j,:));
        microglia_adj_noPLX(i,j)=temp_cor(1,2);
    end
end
%define community based on Louvain method
A_noPLX=abs(microglia_adj_noPLX);
gamma = 1;
k = full(sum(A_noPLX));
twom = sum(k);
B = full(A_noPLX - gamma*k'*k/twom);
[S_noPLX,Q] = genlouvain(B);
com_color_noPLX=zeros(234,3);
for i = 1:234
    com_color_noPLX(i,S_noPLX(i))=1;
end
[S_sorted_noPLX,S_idx_noPLX]=sort(S_noPLX);
nodeColor_noPLX=zeros(length(S_noPLX),3);
for i = 1:length(S_noPLX)-1
    nodeColor_noPLX(i,S_sorted_noPLX(i))=1;
end

norm_microglia_rearrange_noPLX=norm_microglia_gene(S_idx_noPLX,groups_noPLX);
microglia_adj_rearrange_noPLX=[];
for i = 1:size(norm_microglia_rearrange_noPLX,1)
    for j = 1:size(norm_microglia_rearrange_noPLX,1)
        temp_cor=corrcoef(norm_microglia_rearrange_noPLX(i,:),norm_microglia_rearrange_noPLX(j,:));
        microglia_adj_rearrange_noPLX(i,j)=temp_cor(1,2);
    end
end
microglia_adj_rearrange_noPLX = abs(microglia_adj_rearrange_noPLX);
[~,rank_noPLX]=sort(sum(microglia_adj_rearrange_noPLX));
name_microglia_rearrange_noPLX=microglia_name(S_idx_noPLX);
H=schemaball(microglia_adj_rearrange_noPLX,name_microglia_rearrange_noPLX,[1,1,1],[]);
H.s.CData=nodeColor_noPLX(transpose(rank_noPLX),:);
H.s.MarkerEdgeAlpha=1;
H.s.MarkerEdgeColor=[0.1,0.1,0.1];
H.s.LineWidth=0.5;
H.s.MarkerFaceAlpha=0.5;
H.s.SizeData=500.^0.5;

%color the old schemaball with the noPLX group
nodeColor_crossgroup=zeros(length(S_noPLX),3);
colororder=S_idx(rank);
for i=2:length(colororder)
    temp_pos=find(S_idx_noPLX==colororder(i));
    nodeColor_crossgroup(i,S_sorted_noPLX(temp_pos))=1;
end
H=schemaball(microglia_adj_rearrange,name_microglia_rearrange,[1,1,1],[]);
H.s.CData=nodeColor_crossgroup;
H.s.MarkerEdgeAlpha=1;
H.s.MarkerEdgeColor=[0.1,0.1,0.1];
H.s.LineWidth=0.5;
H.s.MarkerFaceAlpha=0.5;
H.s.SizeData=500.^0.5;

%% Venn diagram, figure B
for i=1:size(norm_microglia_gene,1)
    [~,microglia_p(i,1)]=ttest2(norm_microglia_gene(i,21:30),norm_microglia_gene(i,53:58)); %Forebrain
    [~,microglia_p(i,2)]=ttest2(norm_microglia_gene(i,31:40),norm_microglia_gene(i,59:64)); %Hindbrain
end
threshold=0.01;
microglia_sig=microglia_p<threshold;
group1=S==1;
group2=S==2;
group3=S==3;
%total genes
sum(group1)
%modulated genes
sum((microglia_sig(:,1)|microglia_sig(:,2))&group1)
%forebrain modulated genes
sum((microglia_sig(:,1))&group1)
%hindbrain modulated genes
sum((microglia_sig(:,2))&group1)
%intersect
sum((microglia_sig(:,1)&microglia_sig(:,2))&group1)


%total genes
sum(group2)
%modulated genes
sum((microglia_sig(:,1)|microglia_sig(:,2))&group2)
%forebrain modulated genes
sum((microglia_sig(:,1))&group2)
%hindbrain modulated genes
sum((microglia_sig(:,2))&group2)
%intersect
sum((microglia_sig(:,1)&microglia_sig(:,2))&group2)

%total genes
sum(group3)
%modulated genes
sum((microglia_sig(:,1)|microglia_sig(:,2))&group3)
%forebrain modulated genes
sum((microglia_sig(:,1))&group3)
%hindbrain modulated genes
sum((microglia_sig(:,2))&group3)
%intersect
sum((microglia_sig(:,1)&microglia_sig(:,2))&group3)

%% DAM genes panel D
DAM_idx=0;
for i=1:length(DAM_names)
    for j=1:length(microglia_name)
        if matches(string(microglia_name(j)),string(DAM_names(i)))
            DAM_idx(i)=j;
            break
        end
    end
end

DAM_s=S(DAM_idx);
DAM_idx_active=DAM_idx(9:19);
DAM_names_active=DAM_names(9:19);
[~,DAM_idx_tem]=sort(DAM_s(9:19));
DAM_idx_rearrage=[DAM_idx(1:8),DAM_idx_active(DAM_idx_tem)];
DAM_names_rearrage=[DAM_names(1:8);DAM_names_active(DAM_idx_tem)];

DAM_expression_raw=norm_microglia_gene(DAM_idx_rearrage,:);
DAM_expression_ave(:,1)=mean(DAM_expression_raw(:,53:58),2);
DAM_expression_ave(:,2)=mean(DAM_expression_raw(:,21:30),2);
DAM_expression_ave(:,3)=mean(DAM_expression_raw(:,1:10),2);
DAM_expression_ave(:,4)=mean(DAM_expression_raw(:,59:64),2);
DAM_expression_ave(:,5)=mean(DAM_expression_raw(:,31:40),2);
DAM_expression_ave(:,6)=mean(DAM_expression_raw(:,11:20),2);

%% correlation panel E
nonmicroglia_gene=all_raw(~logical(celltype(:,2)),:);
norm_nonmicroglia_gene=nonmicroglia_gene./(max(nonmicroglia_gene,[],2));
fore_group_idx=[1:10,21:30,41:46,53:58];
hind_group_idx=[11:20,31:40,47:52,59:64];
for i=1:length(fore_group_idx)
    for j=1:length(fore_group_idx)
        tempA=norm_nonmicroglia_gene(:,fore_group_idx(i));
        tempB=norm_nonmicroglia_gene(:,fore_group_idx(j));
        temp_corr=corrcoef(tempA,tempB);
        fore_corr(i,j)=temp_corr(1,2);
    end
end

for i=1:length(hind_group_idx)
    for j=1:length(hind_group_idx)
        tempA=norm_nonmicroglia_gene(:,hind_group_idx(i));
        tempB=norm_nonmicroglia_gene(:,hind_group_idx(j));
        temp_corr=corrcoef(tempA,tempB);
        hind_corr(i,j)=temp_corr(1,2);
    end
end

fore_co_PLXtoWT=reshape(fore_corr(1:10,27:32),[],1);
fore_co_VEHtoWT=reshape(fore_corr(11:20,27:32),[],1);
hind_co_PLXtoWT=reshape(hind_corr(1:10,27:32),[],1);
hind_co_VEHtoWT=reshape(hind_corr(11:20,27:32),[],1);

%% PLS decoding and projection
%start with neuron
neuron_gene=all_raw(~logical(celltype(:,2)),:);
norm_neuron_gene=neuron_gene./(max(neuron_gene,[],2));
%fore brain
Xall=tsne(transpose(norm_neuron_gene(:,fore_group_idx)));
Yall=transpose([1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4]);
figure
gscatter(Xall(:,1),Xall(:,2),Yall);

X_raw=norm_neuron_gene(setdiff([1:1599],IEG_idx),:);
X=norm_neuron_gene(:,hind_group_idx([11:20,27:32]));
Y=[1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2];
X_input=transpose(X);
Y_input=transpose(Y);
[XL,~,XS,~,~,PCTVAR] = plsregress(X_input,Y_input,5);
figure
gscatter(XS(:,1),XS(:,2),Y);
plot(1:5,cumsum(100*PCTVAR(2,:)),'-bo');
X0=mean(X_input);
X0_wtplx=[];
X0_tgplx=[];
for i=1:length(X0)
    X0_wtplx(:,i)=transpose(norm_neuron_gene(i,hind_group_idx(21:26)))-X0(i);
end
Xpro_wtplx=(XL\X0_wtplx')';
hold on
scatter(Xpro_wtplx(:,1),Xpro_wtplx(:,2))
hold off
for i=1:length(X0)
    X0_tgplx(:,i)=transpose(norm_neuron_gene(i,hind_group_idx(1:10)))-X0(i);
end
Xpro_tgplx=(XL\X0_tgplx')';
hold on
scatter(Xpro_tgplx(:,1),Xpro_tgplx(:,2))
hold off


PLS_Pro_combine=0;
PLS_Pro_combine=[XS;Xpro_wtplx;Xpro_tgplx];
WT_Veh_vector=mean(XS(11:16,:));
angle=0;
for i=1:size(PLS_Pro_combine,1)
    u=WT_Veh_vector;
    v=PLS_Pro_combine(i,:);
    CosTheta = max(min(dot(u,v)/(norm(u)*norm(v)),1),-1);
    angle(i) = real(acosd(CosTheta));
end
figure
plot(angle)
distance=0;
for i=1:size(PLS_Pro_combine,1)
    distance(i) = mahal(PLS_Pro_combine(i,:),XS(11:16,:));
end
PLS_Pro_combine=[PLS_Pro_combine,transpose([1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4])];
figure
plot(distance)

%hind brain
Xall=tsne(transpose(norm_neuron_gene(:,hind_group_idx)));
Yall=transpose([1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4]);
figure
gscatter(Xall(:,1),Xall(:,2),Yall);

X=norm_neuron_gene(:,hind_group_idx([11:20,27:32]));
Y=[1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2];
X_input=transpose(X);
Y_input=transpose(Y);
[XL,~,XS] = plsregress(X_input,Y_input,5);
figure
gscatter(XS(:,1),XS(:,2),Y);
X0=mean(X_input);
X0_wtplx=[];
X0_tgplx=[];
for i=1:length(X0)
    X0_wtplx(:,i)=transpose(norm_neuron_gene(i,hind_group_idx(21:26)))-X0(i);
end
Xpro_wtplx=(XL\X0_wtplx')';
hold on
scatter(Xpro_wtplx(:,1),Xpro_wtplx(:,2))
hold off
for i=1:length(X0)
    X0_tgplx(:,i)=transpose(norm_neuron_gene(i,hind_group_idx(1:10)))-X0(i);
end
Xpro_tgplx=(XL\X0_tgplx')';
hold on
scatter(Xpro_tgplx(:,1),Xpro_tgplx(:,2))
hold off


PLS_Pro_combine=0;
PLS_Pro_combine=[XS;Xpro_wtplx;Xpro_tgplx];
WT_Veh_vector=mean(XS(11:16,:));
angle=0;
for i=1:size(PLS_Pro_combine,1)
    u=WT_Veh_vector;
    v=PLS_Pro_combine(i,:);
    CosTheta = max(min(dot(u,v)/(norm(u)*norm(v)),1),-1);
    angle(i) = real(acosd(CosTheta));
end
distance=0;
for i=1:size(PLS_Pro_combine,1)
    distance(i) = mahal(PLS_Pro_combine(i,:),XS(11:16,:));
end
PLS_Pro_combine=[PLS_Pro_combine,transpose([1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4])];


%Astrocyte
astro_gene=all_raw(logical(celltype(:,3)),:);
norm_astro_gene=astro_gene./(max(astro_gene,[],2));
%fore brain
Xall=tsne(transpose(norm_astro_gene(:,fore_group_idx)));
Yall=transpose([1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4]);
figure
gscatter(Xall(:,1),Xall(:,2),Yall);

X=norm_astro_gene(:,fore_group_idx([11:20,27:32]));
Y=[1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2];
X_input=transpose(X);
Y_input=transpose(Y);
[XL,~,XS] = plsregress(X_input,Y_input,5);
figure
gscatter(XS(:,1),XS(:,2),Y);
X0=mean(X_input);
X0_wtplx=[];
for i=1:length(X0)
    X0_wtplx(:,i)=transpose(norm_astro_gene(i,fore_group_idx(21:26)))-X0(i);
end
Xpro_wtplx=(XL\X0_wtplx')';
hold on
scatter(Xpro_wtplx(:,1),Xpro_wtplx(:,2))
hold off
X0_tgplx=[];
for i=1:length(X0)
    X0_tgplx(:,i)=transpose(norm_astro_gene(i,fore_group_idx(1:10)))-X0(i);
end
Xpro_tgplx=(XL\X0_tgplx')';
hold on
scatter(Xpro_tgplx(:,1),Xpro_tgplx(:,2))
hold off


PLS_Pro_combine=0;
PLS_Pro_combine=[XS;Xpro_wtplx;Xpro_tgplx];
WT_Veh_vector=mean(XS(11:16,:));
angle=0;
for i=1:size(PLS_Pro_combine,1)
    u=WT_Veh_vector;
    v=PLS_Pro_combine(i,:);
    CosTheta = max(min(dot(u,v)/(norm(u)*norm(v)),1),-1);
    angle(i) = real(acosd(CosTheta));
end
figure
plot(angle)
distance=0;
for i=1:size(PLS_Pro_combine,1)
    distance(i) = mahal(PLS_Pro_combine(i,:),XS(11:16,:));
end
PLS_Pro_combine=[PLS_Pro_combine,transpose([1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4])];
figure
plot(distance)


%hind brain
Xall=tsne(transpose(norm_astro_gene(:,hind_group_idx)));
Yall=transpose([1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4]);
figure
gscatter(Xall(:,1),Xall(:,2),Yall);

X=norm_astro_gene(:,hind_group_idx([11:20,27:32]));
Y=[1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2];
X_input=transpose(X);
Y_input=transpose(Y);
[XL,~,XS] = plsregress(X_input,Y_input,5);
figure
gscatter(XS(:,1),XS(:,2),Y);
X0=mean(X_input);
X0_wtplx=[];
for i=1:length(X0)
    X0_wtplx(:,i)=transpose(norm_astro_gene(i,hind_group_idx(21:26)))-X0(i);
end
Xpro_wtplx=(XL\X0_wtplx')';
hold on
scatter(Xpro_wtplx(:,1),Xpro_wtplx(:,2))
hold off
X0_tgplx=[];
for i=1:length(X0)
    X0_tgplx(:,i)=transpose(norm_astro_gene(i,hind_group_idx(1:10)))-X0(i);
end
Xpro_tgplx=(XL\X0_tgplx')';
hold on
scatter(Xpro_tgplx(:,1),Xpro_tgplx(:,2))
hold off


PLS_Pro_combine=0;
PLS_Pro_combine=[XS;Xpro_wtplx;Xpro_tgplx];
WT_Veh_vector=mean(XS(11:16,:));
angle=0;
for i=1:size(PLS_Pro_combine,1)
    u=WT_Veh_vector;
    v=PLS_Pro_combine(i,:);
    CosTheta = max(min(dot(u,v)/(norm(u)*norm(v)),1),-1);
    angle(i) = real(acosd(CosTheta));
end
figure
plot(angle)
distance=0;
for i=1:size(PLS_Pro_combine,1)
    distance(i) = mahal(PLS_Pro_combine(i,:),XS(11:16,:));
end
PLS_Pro_combine=[PLS_Pro_combine,transpose([1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4])];
figure
plot(distance)

%% PLX upregulation genes
PLX_up=zeros(size(norm_microglia_gene,1),1);
for i=1:size(norm_microglia_gene,1)
    temp_TG_F_PLX=norm_microglia_rearrange(i,1:10);
    temp_TG_F_Veh=norm_microglia_rearrange(i,21:30);
    temp_TG_H_PLX=norm_microglia_rearrange(i,11:20);
    temp_TG_H_Veh=norm_microglia_rearrange(i,31:40);
    if (mean(temp_TG_F_PLX)>mean(temp_TG_F_Veh))&(mean(temp_TG_H_PLX)>mean(temp_TG_H_Veh))
        PLX_up(i)=1;
    end
end
PLX_up_idx=find(PLX_up>0);
temp_TG_F_PLX=norm_microglia_rearrange(PLX_up_idx,1:10);
temp_TG_F_Veh=norm_microglia_rearrange(PLX_up_idx,21:30);
for i=1:length(PLX_up_idx)
    tempA=temp_TG_F_PLX(i,:);
    tempB=temp_TG_F_Veh(i,:);
end
x=[1:length(PLX_up_idx)];
data=[mean(temp_TG_F_PLX,2),mean(temp_TG_F_Veh,2)];
sd=[std(temp_TG_F_PLX,0,2)./sqrt(10),std(temp_TG_F_Veh,0,2)./sqrt(10)];
plx_up_names=name_microglia_rearrange(PLX_up_idx);
x_name=categorical(plx_up_names);
figure
bar(x,data);

hold on
er1 = errorbar(x-0.15,data(:,1),sd(:,1));
er2 = errorbar(x+0.15,data(:,2),sd(:,2));  
er1.Color = [0 0 0];                            
er1.LineStyle = 'none';
er2.Color = [0 0 0];                            
er2.LineStyle = 'none'; 
hold off

temp_TG_H_PLX=norm_microglia_rearrange(PLX_up_idx,11:20);
temp_TG_H_Veh=norm_microglia_rearrange(PLX_up_idx,31:40);
sd=[std(temp_TG_H_PLX,0,2)./sqrt(10),std(temp_TG_H_Veh,0,2)./sqrt(10)];
x=[1:length(PLX_up_idx)];
data=[mean(temp_TG_H_PLX,2),mean(temp_TG_H_Veh,2)];
figure
bar(x,data);
hold on
er1 = errorbar(x-0.15,data(:,1),sd(:,1));
er2 = errorbar(x+0.15,data(:,2),sd(:,2));  
er1.Color = [0 0 0];                            
er1.LineStyle = 'none';
er2.Color = [0 0 0];                            
er2.LineStyle = 'none'; 
hold off

%% male / female difference
neuron_gene=all_raw(~logical(celltype(:,2)),:);
norm_neuron_gene=neuron_gene./(max(neuron_gene,[],2));
%fore brain
X=norm_neuron_gene(:,[21:30,53:58]);
Y=[1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2];
X_input=transpose(X);
Y_input=transpose(Y);
[XL,~,XS,~,~,PCTVAR] = plsregress(X_input,Y_input,5);
%figure
%plot(1:5,cumsum(100*PCTVAR(2,:)),'-bo');
figure('Position',[200,500,900,400]);
subplot(1,8,1:3)
gscatter(XS(:,1),XS(:,2),Y);
X0=mean(X_input);
X0_plxm=[];
X0_plxf=[];
for i=1:length(X0)
    X0_plxm(:,i)=transpose(norm_neuron_gene(i,male_fore_PLX))-X0(i);
end
Xpro_plxm=(XL\X0_plxm')';
hold on
scatter(Xpro_plxm(:,1),Xpro_plxm(:,2))
hold off
for i=1:length(X0)
    X0_plxf(:,i)=transpose(norm_neuron_gene(i,female_fore_PLX))-X0(i);
end
Xpro_plxf=(XL\X0_plxf')';
hold on
scatter(Xpro_plxf(:,1),Xpro_plxf(:,2))
hold off
legend('TG','WT','M+PLX','F+PLX');
xlabel('PLS dimention 1');
ylabel('PLS dimention 2');

PLS_Pro_combine=[Xpro_plxm;Xpro_plxf];
distance=0;
for i=1:size(PLS_Pro_combine,1)
    distance(i) = mahal(PLS_Pro_combine(i,:),XS(11:16,:));
end
subplot(1,8,4)
bar([mean(distance(1:5)),mean(distance(6:10))]);
hold on

er = errorbar([1,2],[mean(distance(1:5)),mean(distance(6:10))],[],[std(distance(1:5)),std(distance(6:10))]);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
hold off
set(gca,'xticklabel',{'M','F'});
ylabel('Distance to WT');

X=norm_neuron_gene(:,[male_fore_PLX,female_fore_PLX,53:58,21:30]);
Y=[1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4];
X_input=transpose(X);
Y_input=transpose(Y);
[XL,~,XS,~,~,PCTVAR] = plsregress(X_input,Y_input,5);
subplot(1,8,6:8)
gscatter(XS(:,1),XS(:,2),Y);
legend('M+PLX','F+PLX','WT','TG');
xlabel('PLS dimention 1');
ylabel('PLS dimention 2');

%hind brain
neuron_gene=all_raw(logical(celltype(:,1)),:);
norm_neuron_gene=neuron_gene./(max(neuron_gene,[],2));
X=norm_neuron_gene(:,[31:40,59:64]);
Y=[1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2];
X_input=transpose(X);
Y_input=transpose(Y);
[XL,~,XS,~,~,PCTVAR] = plsregress(X_input,Y_input,5);
%figure
%plot(1:5,cumsum(100*PCTVAR(2,:)),'-bo');
figure('Position',[200,500,900,400]);
subplot(1,8,1:3)
gscatter(XS(:,1),XS(:,2),Y);
X0=mean(X_input);
X0_plxm=[];
X0_plxf=[];
for i=1:length(X0)
    X0_plxm(:,i)=transpose(norm_neuron_gene(i,male_hind_PLX))-X0(i);
end
Xpro_plxm=(XL\X0_plxm')';
hold on
scatter(Xpro_plxm(:,1),Xpro_plxm(:,2))
hold off
for i=1:length(X0)
    X0_plxf(:,i)=transpose(norm_neuron_gene(i,female_hind_PLX))-X0(i);
end
Xpro_plxf=(XL\X0_plxf')';
hold on
scatter(Xpro_plxf(:,1),Xpro_plxf(:,2))
hold off
legend('TG','WT','M+PLX','F+PLX');
xlabel('PLS dimention 1');
ylabel('PLS dimention 2');

PLS_Pro_combine=[Xpro_plxm;Xpro_plxf];
distance=0;
for i=1:size(PLS_Pro_combine,1)
    distance(i) = mahal(PLS_Pro_combine(i,:),XS(11:16,:));
end
subplot(1,8,4)
bar([mean(distance(1:5)),mean(distance(6:10))]);
hold on

er = errorbar([1,2],[mean(distance(1:5)),mean(distance(6:10))],[],[std(distance(1:5)),std(distance(6:10))]);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
hold off
set(gca,'xticklabel',{'M','F'});
ylabel('Distance to WT');

X=norm_neuron_gene(:,[male_fore_PLX,female_fore_PLX,59:64,31:40]);
Y=[1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4];
X_input=transpose(X);
Y_input=transpose(Y);
[XL,~,XS,~,~,PCTVAR] = plsregress(X_input,Y_input,5);
subplot(1,8,6:8)
gscatter(XS(:,1),XS(:,2),Y);
legend('M+PLX','F+PLX','WT','TG');
xlabel('PLS dimention 1');
ylabel('PLS dimention 2');

%% sex difference single gene level
% Let's do a volcano plot
% fore brain
for i=1:size(all_raw,1)
    tempM = all_raw(i,male_fore_PLX);
    tempF = all_raw(i,female_fore_PLX);
    ratio(i)=mean(tempM)/mean(tempF);
    [~,pvalue(i)]=ttest2(tempM,tempF);
end
scatter(ratio,pvalue);
cand=[];
for i=1:size(all_raw,1)
    if pvalue(i)<0.01 | ratio(i)>1.5
        text(ratio(i),pvalue(i),names_sort{i})
        cand=[cand,i];
    end
end

%hind brain
for i=1:size(all_raw,1)
    tempM = all_raw(i,male_hind_PLX);
    tempF = all_raw(i,female_hind_PLX);
    ratio(i)=mean(tempM)/mean(tempF);
    [~,pvalue(i)]=ttest2(tempM,tempF);
end
scatter(ratio,pvalue);
cand=[];
for i=1:size(all_raw,1)
    if pvalue(i)<0.01 | ratio(i)>1.5
        text(ratio(i),pvalue(i),names_sort{i})
        cand=[cand,i];
    end
end




for i=1:size(all_raw,1)
    tempM = all_raw(i,21:30);
    tempF = all_raw(i,53:58);
    ratio(i)=mean(tempM)/mean(tempF);
    [~,pvalue(i)]=ttest2(tempM,tempF);
end
scatter(ratio,pvalue);
cand=[];
for i=1:size(all_raw,1)
    if pvalue(i)<0.01 | ratio(i)>1.5
        text(ratio(i),pvalue(i),names_sort{i})
        cand=[cand,i];
    end
end

for i=1:size(all_raw,1)
    tempM = all_raw(i,31:40);
    tempF = all_raw(i,59:64);
    ratio(i)=mean(tempM)/mean(tempF);
    [~,pvalue(i)]=ttest2(tempM,tempF);
end
scatter(ratio,pvalue);
cand=[];
for i=1:size(all_raw,1)
    if pvalue(i)<0.01 | ratio(i)>1.5
        text(ratio(i),pvalue(i),names_sort{i})
        cand=[cand,i];
    end
end
%% male female separate in veh group
neuron_gene=all_raw(~logical(celltype(:,2)),:);
norm_neuron_gene=neuron_gene./(max(neuron_gene,[],2));
%fore brain
X=norm_neuron_gene(:,[21:30,53:58]);
Y=[1,2,2,2,2,1,1,2,1,1,2,1,1,2,1,1];
Yall=[1,2,2,2,2,1,1,2,1,1,4,3,3,4,3,3];
X_input=transpose(X);
Y_input=transpose(Y);
Xall=tsne(X_input);
gscatter(Xall(:,1),Xall(:,2),Y);
[XL,~,XS,~,~,PCTVAR] = plsregress(X_input,Y_input,5);
figure
%plot(1:5,cumsum(100*PCTVAR(2,:)),'-bo');
%figure('Position',[200,500,900,400]);
%subplot(1,8,1:3)
gscatter(XS(:,1),XS(:,2),Yall);
X0=mean(X_input);
X0_plxm=[];
X0_plxf=[];
for i=1:length(X0)
    X0_plxm(:,i)=transpose(norm_neuron_gene(i,male_fore_PLX))-X0(i);
end
Xpro_plxm=(XL\X0_plxm')';
hold on
scatter(Xpro_plxm(:,1),Xpro_plxm(:,2))
hold off
for i=1:length(X0)
    X0_plxf(:,i)=transpose(norm_neuron_gene(i,female_fore_PLX))-X0(i);
end
Xpro_plxf=(XL\X0_plxf')';
hold on
scatter(Xpro_plxf(:,1),Xpro_plxf(:,2))
hold off
legend('TG_M','TG_F','WT_M','WT_F','M+PLX','F+PLX');
xlabel('PLS dimention 1');
ylabel('PLS dimention 2');

figure
plot(XL(:,2));
names_nonmicroglia=names_sort(~logical(celltype(:,2)),:);
for i=1:size(XL,1)
    if abs(XL(i,2))>0.4
        text(i,XL(i,2),names_nonmicroglia{i})
    end
end





%hindbrain
X=norm_neuron_gene(:,[31:40,59:64]);
Y=[1,2,2,2,2,1,1,2,1,1,2,1,1,2,1,1];
Yall=[1,2,2,2,2,1,1,2,1,1,4,3,3,4,3,3];
X_input=transpose(X);
Y_input=transpose(Y);
Xall=tsne(X_input);
gscatter(Xall(:,1),Xall(:,2),Y);
[XL,~,XS,~,~,PCTVAR] = plsregress(X_input,Y_input,5);
figure
%plot(1:5,cumsum(100*PCTVAR(2,:)),'-bo');
gscatter(XS(:,1),XS(:,2),Yall);
X0=mean(X_input);
X0_plxm=[];
X0_plxf=[];
for i=1:length(X0)
    X0_plxm(:,i)=transpose(norm_neuron_gene(i,male_fore_PLX))-X0(i);
end
Xpro_plxm=(XL\X0_plxm')';
hold on
scatter(Xpro_plxm(:,1),Xpro_plxm(:,2))
hold off
for i=1:length(X0)
    X0_plxf(:,i)=transpose(norm_neuron_gene(i,female_fore_PLX))-X0(i);
end
Xpro_plxf=(XL\X0_plxf')';
hold on
scatter(Xpro_plxf(:,1),Xpro_plxf(:,2))
hold off
legend('TG_M','TG_F','WT_M','WT_F','M+PLX','F+PLX');
xlabel('PLS dimention 1');
ylabel('PLS dimention 2');

figure
plot(XL(:,2));
names_nonmicroglia=names_sort(~logical(celltype(:,2)),:);
for i=1:size(XL,1)
    if abs(XL(i,2))>0.4
        text(i,XL(i,2),names_nonmicroglia{i})
    end
end


% 3D TG vs WT ;  M vs F; Hind vs Fore;
X=norm_neuron_gene(:,[21:40,53:64,71:76]);
sex=[1,2,2,2,2,1,1,2,1,1,1,2,2,2,2,1,1,2,1,1,2,1,1,2,1,1,2,1,1,2,1,1,1,2,2,1,2,2];
transgene=[ones(1,20),zeros(1,18)];
region=[1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,1,1,1,1,1,1,2,2,2,2,2,2,1,1,1,2,2,2];
Y=[sex;transgene;region];
Yall=[1,2,2,2,2,1,1,2,1,1,5,6,6,6,6,5,5,6,5,5,4,3,3,4,3,3,8,7,7,8,7,7,3,4,4,7,8,8];
X_input=transpose(X);
Y_input=transpose(Y);
[XL,~,XS,~,~,PCTVAR] = plsregress(X_input,Y_input,5);
figure
%plot(1:5,cumsum(100*PCTVAR(2,:)),'-bo');
gscatter3(XS(:,1),XS(:,2),XS(:,3),Yall);


% 3D PLX conc, M vs F
% fore brain only, tg only
neuron_gene=all_raw(~logical(celltype(:,2)),:);
norm_neuron_gene=neuron_gene./(max(neuron_gene,[],2));
X=norm_neuron_gene(:,[11:20,31:40]);
X=norm_neuron_gene(IEG_idx,[11:20,31:40]);
sex=[1,1,2,1,1,2,2,2,1,2,1,2,2,2,2,1,1,2,1,1];
PLX=transpose(PLX_conc([11:20,31:40]));
Yall=[1,1,2,1,1,2,2,2,1,2,3,4,4,4,4,3,3,4,3,3];
Y=[sex;PLX];
X_input=transpose(X);
Y_input=transpose(Y);
[XL,yl,XS,YS,beta,PCTVAR,MSE,stats] = plsregress(X_input,Y_input,5);
W0 = stats.W ./ sqrt(sum(stats.W.^2,1));
p = size(XL,1);
sumSq = sum(XS.^2,1).*sum(yl.^2,1);
vipScore = sqrt(p* sum(sumSq.*(W0.^2),2) ./ sum(sumSq,2));
vip_idx=find(vipScore>1.5);
plot(abs(W0(vip_idx,1:2)));

plotsubgroup(neuron_gene,vip_idx,neuron_name);

figure
%plot(1:5,cumsum(100*PCTVAR(2,:)),'-bo');
gscatter3(XS(:,1),XS(:,2),XS(:,3),Yall);