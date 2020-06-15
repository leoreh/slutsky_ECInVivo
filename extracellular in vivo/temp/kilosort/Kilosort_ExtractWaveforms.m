function waveforms_all = Kilosort_ExtractWaveforms(rez)
        % Extracts waveforms from a dat file using GPU enable filters.
        % Based on the GPU enable filter from Kilosort.
        % All settings and content are extracted from the rez input structure
        %
        % Inputs:
        %   rez -  rez structure from Kilosort
        %
        % Outputs:
        %   waveforms_all - structure with extracted waveforms
        
        % Extracting content from the .rez file
        ops = rez.ops;
        NT = ops.NT;
        if exist('ops.fbinary') == 0
            warning(['Binary file does not exist: ', ops.fbinary])
        end
        d = dir(ops.fbinary);

        NchanTOT = ops.NchanTOT;
        chanMap = ops.chanMap;
        chanMapConn = chanMap(rez.connected>1e-6);
        kcoords = ops.kcoords;
        ia = rez.ia;
        spktimes = rez.st3(:,1);
        
        if ispc
            dmem         = memory;
            memfree      = dmem.MemAvailableAllArrays/8;
            memallocated = min(ops.ForceMaxRAMforDat, dmem.MemAvailableAllArrays) - memfree;
            memallocated = max(0, memallocated);
        else
            memallocated = ops.ForceMaxRAMforDat;
        end
        ops.ForceMaxRAMforDat   = 10000000000;
        memallocated = ops.ForceMaxRAMforDat;
        nint16s      = memallocated/2;
        
        NTbuff      = NT + 4*ops.ntbuff;
        Nbatch      = ceil(d.bytes/2/NchanTOT /(NT-ops.ntbuff));
        Nbatch_buff = floor(4/5 * nint16s/ops.Nchan /(NT-ops.ntbuff)); % factor of 4/5 for storing PCs of spikes
        Nbatch_buff = min(Nbatch_buff, Nbatch);
        
        DATA =zeros(NT, NchanTOT,Nbatch_buff,'int16');
        
        if isfield(ops,'fslow')&&ops.fslow<ops.fs/2
            [b1, a1] = butter(3, [ops.fshigh/ops.fs,ops.fslow/ops.fs]*2, 'bandpass');
        else
            [b1, a1] = butter(3, ops.fshigh/ops.fs*2, 'high');
        end
        
        if isfield(ops,'xml')
            disp('Loading xml from rez for probe layout')
            xml = ops.xml;
        elseif exist(fullfile(ops.root,[ops.basename,'.xml']))==2
            disp('Loading xml for probe layout from root folder')
            xml = LoadXml(fullfile(ops.root,[ops.basename,'.xml']));
            ops.xml = xml;
        end
        
        fid = fopen(ops.fbinary, 'r');
        
        waveforms_all = [];
%         kcoords2 = unique(ops.kcoords);
        unitGrp = kcoords(ampMaxCh);
        kcoords2 = unique(unitGrp);

        channel_order = {};
        indicesTokeep = {};
%         connected_index = zeros(size(rez.connected));
%         connected_index(rez.connected)=1:length(chanMapConn);
        
        for i = 1:length(kcoords2)
            grp = kcoords2(i);
            waveforms_all{i} = zeros(sum(kcoords==grp),ops.nt0,size(rez.ia{i},1));
            if exist('xml')
                [channel_order,channel_index] = sort(xml.AnatGrps(kcoords2(i)).Channels+1);
                [~,indicesTokeep{i},~] = intersect(chanMapConn,channel_order);
                
                %indicesTokeep{i} = connected_index(indicesTokeep{i});
            end
        end
        
        fprintf('Extraction of waveforms begun \n')
        for ibatch = 1:Nbatch
            if mod(ibatch,10)==0
                if ibatch~=10
                    fprintf(repmat('\b',[1 length([num2str(round(100*(ibatch-10)/Nbatch)), ' percent complete'])]))
                end
                fprintf('%d percent complete', round(100*ibatch/Nbatch));
            end
            
            offset = max(0, 2*NchanTOT*((NT - ops.ntbuff) * (ibatch-1) - 2*ops.ntbuff));
            if ibatch==1
                ioffset = 0;
            else
                ioffset = ops.ntbuff;
            end
            fseek(fid, offset, 'bof');
            buff = fread(fid, [NchanTOT NTbuff], '*int16');
            
            %         keyboard;
            
            if isempty(buff)
                break;
            end
            nsampcurr = size(buff,2);
            if nsampcurr<NTbuff
                buff(:, nsampcurr+1:NTbuff) = repmat(buff(:,nsampcurr), 1, NTbuff-nsampcurr);
            end
            if ops.GPU
                dataRAW = gpuArray(buff);
            else
                dataRAW = buff;
            end
            
            dataRAW = dataRAW';
            dataRAW = single(dataRAW);
            dataRAW = dataRAW(:, chanMapConn);
            dataRAW = dataRAW-median(dataRAW,2);
            datr = filter(b1, a1, dataRAW);
            datr = flipud(datr);
            datr = filter(b1, a1, datr);
            datr = flipud(datr);
            DATA = gather_try(int16( datr(ioffset + (1:NT),:)));
            dat_offset = offset/NchanTOT/2+ioffset;
            % Saves the waveforms occuring within each batch
            for i = 1:length(kcoords2)
                grp = kcoords2(i);
%                 ch_subset = 1:length(chanMapConn);
                temp = find(ismember(spktimes(ia{i}), round([ops.nt0/2+1:size(DATA,1)-ops.nt0/2]) + dat_offset));
                temp2 = spktimes(ia{i}(temp))-dat_offset;
                
                startIndicies = temp2-ops.nt0/2+1;
                stopIndicies = temp2+ops.nt0/2;
                X = cumsum(accumarray(cumsum([1;stopIndicies(:)-startIndicies(:)+1]),[startIndicies(:);0]-[0;stopIndicies(:)]-1)+1);
                X = round(X(1:end-1));
                waveforms_all{i}(:,:,temp) = reshape(DATA(X,indicesTokeep{i})',size(indicesTokeep{i},1),ops.nt0,[]);
            end
        end