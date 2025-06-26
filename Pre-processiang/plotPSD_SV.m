%% Use micro g^2/Hz as the unit of PSD and SV
function plotPSD_SV(ichan_plt,ichan_cal,fres,fb,fc,bandselect)
% tdta: time domain data. [nt,n]， nt: length of data  n：number of channels
% fs: sampling rate.
% fres: resolution ratio (0~1). default: 0.1（ref: bayoma sec:4.6）
% fc: upper limit of fre in the figure.
% fb: lower limit of fre in the figure.
% bandselection: 0: don't selet the band; 1: select the band
if isempty(fres)
    fres = 0.1;
end
if isempty(fb)
    fb = 0;
end
[fnames,paths]=uigetfile({'*.mat'},'Pick data file(s)','MultiSelect','on');
fnames = fullfile(paths,fnames);
in = load(fnames);
tdata = in.tdata;
[nt,n] = size(tdata);
if isempty(ichan_cal)
    ichan_cal = 1:n; % default
end
tdata = tdata(:,ichan_cal);
fs = in.fs;
if isempty(fc)
    fc = fs/2;
end
if isfield(in,'time')
    time = in.time;
else
    time = (1:nt)/fs;
end
if isempty(ichan_plt)
    ichan_plt = [1,2]; % default
end
tdata = detrend(tdata);
nplt = length(ichan_plt);
% plot time history and ask user to select time span
%===========================================================
fprintf('Select 2 points defining time segment for analysis.');
figure(466);
for ii = 1:nplt
    subplot(nplt,1,ii);
    jj = ichan_plt(ii); % the actual dof no.
    plot(time,10^6*tdata(:,jj),'LineWidth',1.1);
    if isfield(in,'time')
        datetickzoom('x','HH:MM');
        xlabel('$\rm\bf{Time}\ \left[\it{hh:mm}\right]$','interpreter','latex');
    end
%     ylabel('$\rm\bf{Wind\ velocity\ \left[\it{\mu{m/s}}\right]}$','interpreter','latex');
    ylabel(['$\rm\bf{DOF}\ ','\rm\bf{',num2str(jj),'}\ ','\left[\it{\mu{g}}\right]$'],'interpreter','latex');
    set(gca,'FontSize',12,'FontSmoothing','on','FontWeight','bold');
    set(0,'defaultfigurecolor','w');
%     ylim([-1000 1000]);
    drawnow;
end
[t1t2, ~] = ginput(2);
if length(t1t2)<2
    I1I2 = [1 nt];
else
    cut1 = find(time <= t1t2(1), 1, 'last' );
    if isempty(cut1)
        cut1=1;
    end
    cut2 = find(time <= t1t2(2), 1, 'last' );
    if isempty(cut2) || cut2<=cut1
        cut2=nt;
    end
    I1I2 = [cut1 cut2];
end

close(figure(466));

tdata = tdata(I1I2(1):I1I2(2),:);
[nt,n] = size(tdata);
fprintf(['Selected duration is ',num2str(nt/fs),'sec.\n']);

if fres < fs/nt
    error(['frequency resolution is no finer than ',num2str(fs/nt)]);
end
navg = nt/(fs/fres);
[sw,w]=specden(tdata,1/fs,navg); % averaged spectra
sw = sw(2:end,:,:); w = w(2:end); % skip D.C.
f = w/2/pi;
nf = length(f);
[~,ind] = min(abs(f-fc));

swtmp = zeros(n,nf);
for ii = 1:n
    swtmp(ii,:) = sw(:,ii,ii);
end
maxpsd = max(max(swtmp(:,1:ind))); minpsd = min(min(swtmp(:,1:ind)));

if n > 5; n = 5; end
Sv2 = zeros(nf,n);
for ii=1:nf
    Sv2(ii,:) = svds(squeeze(real(sw(ii,:,:))),n).';
end
maxSv2 = max(max(Sv2(1:ind,:))); minSv2 = min(min(Sv2(1:ind,:)));

figure;
semilogy(f,10^12*swtmp,'-','LineWidth',1.1);
xlim([fb fc]);
ylim([0.55*10^12*minpsd 10^12*maxpsd*1.5]);
grid on;
xlabel('$\rm\bf{Frequency}\ \left[\it{Hz}\right]$','interpreter','latex');
% ylabel('$\rm\bf{PSD}\ \left[\it{{\left(\mu{m/s}\right)}^{2}/Hz}\right]$','interpreter','latex');
ylabel('$\rm\bf{PSD}\ \left[\it{{\left(\mu{g}\right)}^{2}/Hz}\right]$','interpreter','latex');
set(gca,'FontSize',16,'FontSmoothing','on','FontWeight','bold');
set(0,'defaultfigurecolor','w');

figure;
semilogy(f,10^12*Sv2,'-','LineWidth',1.1);
xlim([fb fc]);
ylim([0.5*10^12*minSv2 1.5*maxSv2*10^12]);
grid on;
xlabel('$\rm\bf{Frequency}\ \left[\it{Hz}\right]$','interpreter','latex');
% ylabel('$\rm\bf{PSD}\ \left[\it{{\left(\mu{m/s}\right)}^{2}/Hz}\right]$','interpreter','latex');
ylabel('$\rm\bf{SV}\ \left[\it{{\left(\mu{g}\right)}^{2}/Hz}\right]$','interpreter','latex');
set(gca,'FontSize',16,'FontSmoothing','on','FontWeight','bold');
set(0,'defaultfigurecolor','w');

hold on
if bandselect==1
    conin1=input('start selecting the band?(1/0) \n');
    if logical(conin1)
        nband=1;
        while nband<100
            [x,~] = ginput;
            n = length(x);
            f1f2(nband,:) = [x(1),x(2)];
            f0{nband,1} = x(3:n)';
            
            yband=ones(1,2)*0.8*10^12*minSv2;
            [rnf,cnf]=size(f0{nband});
            ynf=ones(rnf,cnf)*0.8*10^12*minSv2;
            h1{nband}=plot(f1f2(nband,:)',yband','k','LineWidth',2);
            h2{nband}=plot(f0{nband,1},ynf,'ko','LineWidth',2,'MarkerSize',6);
            h3{nband}=text(f1f2(nband,1),0.8*10^12*minSv2,'[','LineStyle','none','FontWeight','bold');
            h4{nband}=text(f1f2(nband,2),0.8*10^12*minSv2,']','LineStyle','none','FontWeight','bold');
            hold on
            nband=nband+1;
            conin2=input('Please input the option:\n 0: Abandon the band\n 1: The next band\n 2: Reselect the band\n 3: End the selection \n');
            if conin2==0 || conin2==2
                nband=nband-1;
                f1f2(nband,:) = [];
                f0{nband,1} = [];
                delete(h1{nband});delete(h2{nband});delete(h3{nband});delete(h4{nband});
                if conin2==0
                    for nab = nband-1:-1:1
                        conin3=input('Please input the option:\n 0: Abandon the band\n 1: The next band\n 2: Reselect the band\n 3: End the selection \n');
                        if conin3==1
                            break
                        elseif conin3==3
                            conin2=3;
                            break
                        elseif conin3==2
                            f1f2(nab,:) = [];
                            f0{nab,1} = [];
                            delete(h1{nab});delete(h2{nab});delete(h3{nab});delete(h4{nab});
                            break
                        elseif conin3==0
                            f1f2(nab,:) = [];
                            f0{nab,1} = [];
                            delete(h1{nab});delete(h2{nab});delete(h3{nab});delete(h4{nab});
                        end
                    end
                end
            end
            if conin2==3
                break
            end
        end
        save(uiputfile('*.mat'),'f1f2','f0','tdata','fs');
    end
end
end