%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% range of spike counts of sust vs onset

fid=fopen('AC_spkcnt_data.m');
%4313 u513 data is no good, 4671

% filenum, unit, ch, cf, stim_num, fmod_list, 
%    0 stimstart stimstart+50 stimstart+100 stimstart+stimlength stimstart+stimlength+100 stimstart+stimlength+post_rec
% spk_counts organized as nspk(fmod(1), stimstart-0) nspk(fmod(1), stimstart+50-stimstart) ... ... nspk(fmod(2), stimstart-0) nspk(fmod(2), stimstart+50-stimstart)
%

quit=0;
numwin=6; % number of analysis window markers e.g., 0-stimstart, stimstart-stimstart+50, ...
badvalue=nan;
histogrambin=0:.5:15;

pause off

result=[];
result2=[];
cf_result=[];
while ~quit
   s=fgets(fid);
   if s==-1
      quit=1;
      break;
   end
   data=sscanf(s,'%f');
   filenum=data(1);
   unit=data(2);
   ch=data(3);
   cf=data(4);
   numstim=data(5);
   fmodlist=data(6:(6+numstim-1));
   spkcnt=[];
   for i=1:numstim
      spkcnt=[spkcnt; data(5+numstim+(1:numwin)+(6*(i-1)))'];
   end
   
   w1=(spkcnt(:,1)); % 0 - stimstart
   w2=(spkcnt(:,2)); % stimstart - stimstart+50
   w3=(spkcnt(:,3)); % stimstart+50 - stimstart+100
   w4=(spkcnt(:,4)); % stimstart+100 - stimstart+stimlength
   w5=(spkcnt(:,5)); % stimstart+stimlength - stimstart+stimlength+100
   w6=(spkcnt(:,6)); % stimstart+stimlength+100 - stimstart+stimlength+postrec
   
   
   %%%%%% total %%%%%%%%%%%
   rmtf=(w2+w3+w4+w5)/(1.100*10)-w1/(0.500*10);
   fmodlist=data(6:(6+numstim-1));
   tmp=[fmodlist rmtf];
   tmp=sortrows(tmp);
   fmodlist=tmp(:,1);
   rmtf=tmp(:,2);
   [y yi]=max(rmtf);
%   semilogx(fmodlist,rmtf, '-*')
   if y>5
       bmf=fmodlist(yi);
       if bmf>0
          hbmf=y/2;
          hrmtf=rmtf-hbmf; % half bw is at zero
          %find half-bw of low side
          lowval=[];
          highval=[];
          x1=max(find(hrmtf(1:yi)<0));
          if ~isempty(x1)
             lowval=roots(polyfit([fmodlist(x1) fmodlist(x1+1)],[hrmtf(x1) hrmtf(x1+1)],1));
%           else 
%              lowval = nan;
          end

          %find half-bw of high side
          tmphrmtf=hrmtf(yi:end);
          tmpfl=fmodlist(yi:end);
          x1=min(find(tmphrmtf<0));
          if ~isempty(x1) & x1>1
             highval=roots(polyfit([tmpfl(x1) tmpfl(x1-1)],[tmphrmtf(x1) tmphrmtf(x1-1)],1));
             tmpfl(x1);         
    %         line([lowval highval], [hbmf hbmf]);
%           else 
%              highval = nan;
          end	
       end

        if ~isempty(highval-lowval)
            totalbw=[highval-lowval];
            result=[result; filenum unit cf totalbw bmf lowval highval cf/lowval*1000 cf/highval*1000];	
        else
            totalbw=badvalue;
            %bmf=badvalue;
        end		
	else
		totalbw=badvalue;
		bmf=badvalue;
	end		



end	
fclose(fid);

%% Plot Resolvability
figure
hold on
for i = 1:length(result)
    plot(result(i,8:9), result(i,6:7)); 
end
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');


