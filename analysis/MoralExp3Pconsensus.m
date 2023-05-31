%% Acá están las densidades y varianzas por grupo, y el gráfico de proba de consenso vs. densidad.

clearvars
%cd 'E:\Escritorio\Todo\Nuevo Análisis'
%cd 'D:\Doctorado\Nuevo Análisis'
%cd 'E:\Nuevo Análisis'

%Datos Originales
% load('DatosMoralS1ConPalabras')
% load('datosS1OrdenPalabras')
% load('datosS1OrdenPartics')
% load('datosS1OrdenPreguntas.mat')
% 
% c=importdata('consensos_moral4.csv');
% consensos=c.data;
% clear c

%Datos Actualizados:
load('DatosMoralS1ConPalabrasUpdated')
% load('datosS1OrdenPalabrasUpdated')
% load('datosS1OrdenParticsUpdated')
% load('datosS1OrdenPreguntasUpdated')

c=importdata('consensos_moral_updated.csv');
consensos=c.data;
clear c


%Promedio de palabras por número de intervención:
% figure()
% for i=9:14
% subplot(3,2,i-8)
% serie=nanmean(datosS1OrdenPalabras(nuevoOrdenPreguntas==i,:),1);
% plot(serie(1:end-1),'o')
% title(['q=',num2str(i)])
% end
% 
% figure()
% serie=nanmean(datosS1OrdenPalabras,1);
% plot(serie(1:end-1),'o')


%%




%%

clearvars
%cd 'E:\Escritorio\Todo\Nuevo Análisis'
cd 'D:\Doctorado\Nuevo Análisis'
%cd 'E:\Nuevo Análisis'

load('DatosMoralS1ConPalabrasUpdated')
load('DatosMoralS3ConPalabrasUpdated')

% load('datosS1OrdenPalabrasUpdated')
% load('datosS1OrdenParticsUpdated')
% load('datosS1OrdenPreguntasUpdated')

c=importdata('consensos_moral_updated.csv');
consensos=c.data;
clear c

%ConTiempos=importdata('featuresMoralTiempos.csv',',');
featuresMoralTiempos = importfileMoral2('featuresMoralTiempos2Updated.csv', 2, 30460);



% Longitud de Palabras (no temporal)


clearvars -except consensos datosS1ConPalabras featuresMoralTiempos

Consensos={[],[],[],[],[],[]};
Rangos={[],[],[],[],[],[]};
Palabras={[],[],[],[],[],[]};
Longitudes={[],[],[],[],[],[]};
PalabrasFacu={[],[],[],[],[],[]};

for q=1:6
    grupos=consensos(consensos(:,2)==q+8,1);
    for grupo=grupos'
        
        consensosPreguntas=consensos(consensos(:,1)==grupo,2);
        consensosRespuestas=consensos(consensos(:,1)==grupo,3);
        consenso=consensosRespuestas(consensosPreguntas==q+8)~=-1;        
        
        if ~isempty(consenso) && ~isempty(featuresMoralTiempos.n_words(featuresMoralTiempos.group_id==grupo & featuresMoralTiempos.question_id==q+8))            
        
        opiniones=datosS1ConPalabras(datosS1ConPalabras(:,7)==grupo,q);
        rangointerno=max(opiniones)-min(opiniones);        
            
        Consensos{q}=[Consensos{q},consenso];
        Rangos{q}=[Rangos{q},rangointerno];
        
        
        PalabrasTiempos=featuresMoralTiempos.n_words(featuresMoralTiempos.group_id==grupo & featuresMoralTiempos.question_id==q+8);
        
        Longitudes{q}=[Longitudes{q},sum(featuresMoralTiempos.words_length__mean(featuresMoralTiempos.group_id==grupo & featuresMoralTiempos.question_id==q+8).*PalabrasTiempos)];
        PalabrasFacu{q}=[PalabrasFacu{q},sum(featuresMoralTiempos.n_words(featuresMoralTiempos.group_id==grupo & featuresMoralTiempos.question_id==q+8))];
        
        end
    end
end


%

TodasPalabras=[PalabrasFacu{:}];
% TodosTurnos=[Turnos{:}];
TodosRangos=[Rangos{:}];
TodosConsensos=[Consensos{:}];
TodasLongitudes=[Longitudes{:}];

%figure()
palabrasRangoConsenso=TodasPalabras(TodosConsensos==1);
palabrasRangoDisenso=TodasPalabras(TodosConsensos==0);
% TurnosRangoConsenso=TodosTurnos(TodosRangos==rangos & TodosConsensos==1);
% TurnosRangoDisenso=TodosTurnos(TodosRangos==rangos & TodosConsensos==0);
longitudesxPalabraConsenso=TodasLongitudes(TodosConsensos==1)./TodasPalabras(TodosConsensos==1);
longitudesxPalabraDisenso=TodasLongitudes(TodosConsensos==0)./TodasPalabras(TodosConsensos==0);


%

figure();
for caso=1:2
    if caso==1
        rangomax=5;%2:10
        rangomin=0;
        ini=2.5;
        ending=5;        
        paso=0.11;
    elseif caso==2
        rangomax=10;
        rangomin=6;
        ini=2.5;
        ending=5.5;
        paso=0.12;
    end

    longitudesxPalabraConsensoRango=TodasLongitudes(TodosRangos<=rangomax & TodosRangos>=rangomin & TodosConsensos==1)./TodasPalabras(TodosRangos<=rangomax & TodosRangos>=rangomin & TodosConsensos==1);
    longitudesxPalabraDisensoRango=TodasLongitudes(TodosRangos<=rangomax & TodosRangos>=rangomin & TodosConsensos==0)./TodasPalabras(TodosRangos<=rangomax & TodosRangos>=rangomin & TodosConsensos==0);


k=1;
cuentas=[];
for longitud=ini:paso:ending
    gruposConsensoDentro=sum(longitudesxPalabraConsensoRango>=longitud & longitudesxPalabraConsensoRango<longitud+paso);
    gruposDisensoDentro=sum(longitudesxPalabraDisensoRango>=longitud & longitudesxPalabraDisensoRango<longitud+paso);
    cuentas(k)=gruposConsensoDentro./(gruposConsensoDentro+gruposDisensoDentro);k=k+1;
end
cuentas
hold on
if caso==1
    %valores=find(cuentas==1.0);
    %cuentas(valores(end))=nan;
    cuentas
    %datos originales
%     plot(0.5:paso:19-paso*3,cuentas(1:end-2),'o','color',[0, 0.5, 0],'LineWidth',1.5)
%     mdl = fitglm(0.5:paso:19-paso*3,cuentas(1:end-2),'Distribution','binomial');
%     
    rangoploteo=ini:paso:ending-3*paso;
    plot(rangoploteo([1:6,8:length(rangoploteo)]),cuentas([1:6,8:length(cuentas)-3]),'o','color',[0, 0.5, 0],'LineWidth',1.5)
    mdl = fitglm(rangoploteo([1:6,8:length(rangoploteo)]),cuentas([1:6,8:length(cuentas)-3]),'Distribution','binomial');
    params = mdl.Coefficients.Estimate;
    hold on
    plot(ini-0.1:paso:ending+paso,1./(exp(-(params(1)+params(2).*(ini-0.1:paso:ending+paso)))+1),'color',[0, 0.5, 0],'LineWidth',1.5)
elseif caso==2
    plot(ini:paso:ending-4*paso,cuentas(1:end-4),'o','color',[0.6350, 0.0780, 0.1840],'LineWidth',1.5)
    mdl = fitglm(ini:paso:ending-4*paso,cuentas(1:end-4),'Distribution','binomial');
    params = mdl.Coefficients.Estimate;
    hold on
    plot(ini-0.1:paso:ending,1./(exp(-(params(1)+params(2).*(ini-0.1:paso:ending)))+1),'color',[0.6350, 0.0780, 0.1840],'LineWidth',1.5)        
end
%lsline
xlabel('Length per Word')
ylabel('P(Consensus)')
xlim([2.4,5])
ylim([0,1])
xticks(2.5:0.5:5)
yticks(0:0.2:1)
legend('Rangos<=5','ajuste','Rangos>5','ajuste','Location','SouthWest')
end    

