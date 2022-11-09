classdef fidnav_b0_retro_recon < handle
  % Linting warning suppression:
  %#ok<*INUSD>  Input argument '' might be unused.  If this is OK, consider replacing it by ~
  %#ok<*NASGU>  The value assigned to variable '' might be unused.
  %#ok<*INUSL>  Input argument '' might be unused, although a later one is used.  Ronsider replacing it by ~
  %#ok<*AGROW>  The variable '' appear to change in size on every loop  iteration. Consider preallocating for speed.
  
  methods
    function process(obj, connection, config, metadata, logging)
      logging.info('Config: \n%s', config);
      
      % Metadata should be MRD formatted header, but may be a string
      % if it failed conversion earlier
      try
        logging.info("Incoming dataset contains %d encodings", numel(metadata.encoding))
        logging.info("First encoding is of type '%s', with field of view of (%g x %g x %g)mm^3, matrix size of (%g x %g x %g), reconstructed size of (%g x %g x %g) and %g coils", ...
          metadata.encoding(1).trajectory, ...
          metadata.encoding(1).encodedSpace.fieldOfView_mm.x, ...
          metadata.encoding(1).encodedSpace.fieldOfView_mm.y, ...
          metadata.encoding(1).encodedSpace.fieldOfView_mm.z, ...
          metadata.encoding(1).encodedSpace.matrixSize.x, ...
          metadata.encoding(1).encodedSpace.matrixSize.y, ...
          metadata.encoding(1).encodedSpace.matrixSize.z, ...
          metadata.encoding(1).reconSpace.matrixSize.x, ...
          metadata.encoding(1).reconSpace.matrixSize.y, ...
          metadata.encoding(1).reconSpace.matrixSize.z, ...
          metadata.acquisitionSystemInformation.receiverChannels)
      catch
        logging.info("Improperly formatted metadata: \n%s", metadata)
      end
      
      % Continuously parse incoming data parsed from MRD messages
      acqGroup = cell(1,0); % ismrmrd.Acquisition;
      navGroup = cell(1,0);
      
      try
        while true
          item = next(connection);
          
          % ----------------------------------------------------------
          % Raw k-space data messages
          % ----------------------------------------------------------
          if isa(item, 'ismrmrd.Acquisition')
            % Accumulate all imaging readouts in a group
            if (~item.head.flagIsSet(item.head.FLAGS.ACQ_IS_NOISE_MEASUREMENT)    && ...
                ~item.head.flagIsSet(item.head.FLAGS.ACQ_IS_PHASECORR_DATA)       && ...
                ~item.head.flagIsSet(item.head.FLAGS.ACQ_IS_PARALLEL_CALIBRATION) && ...
                ~item.head.flagIsSet(item.head.FLAGS.ACQ_IS_NAVIGATION_DATA) )
              acqGroup{end+1} = item;
            end
           
            % Accumulate all imaging readouts in a group
            if (item.head.flagIsSet(item.head.FLAGS.ACQ_IS_NAVIGATION_DATA) )
              navGroup{end+1} = item;
            end
           
            % When this criteria is met, run process_raw() on the accumulated
            % data, which returns images that are sent back to the client.
           % if item.head.flagIsSet(item.head.FLAGS.ACQ_LAST_IN_SLICE)
           if item.head.flagIsSet(item.head.FLAGS.ACQ_LAST_IN_MEASUREMENT)
              logging.info("Processing a group of %d k-space data", length(acqGroup))
              images = obj.process_raw(acqGroup, navGroup, config, metadata, logging);
              logging.debug("Sending image to client")
              connection.send_image(images);
              acqGroup = {};
            end
            
          elseif isempty(item)
            break;
            
          else
            logging.error("Unhandled data type: %s", class(item))
          end
        end
      catch ME
        logging.error(sprintf('%s\nError in %s (%s) (line %d)', ME.message, ME.stack(1).('name'), ME.stack(1).('file'), ME.stack(1).('line')));
      end
      
      % Process any remaining groups of raw or image data.  This can
      % happen if the trigger condition for these groups are not met.
      % This is also a fallback for handling image data, as the last
      % image in a series is typically not separately flagged.
       if ~isempty(acqGroup)
         logging.info("Processing a group of k-space data (untriggered)")
         image = obj.process_raw(acqGroup, navGroup, config, metadata, logging);
         logging.debug("Sending image to client")
         connection.send_image(image);
         acqGroup = cell(1,0);
       end
      
      connection.send_close();
      return
    end
   
    
    function images = process_raw(obj, group, navgroup, config, metadata, logging)
      
        % load calibration data
        load('/tmp/calibrationMatrix.mat'); 
        params_cal = params;
        
        permutedims = [1 2 3];

        swapdims = [0 1]; % ap
      
        params = struct();

        params.NSlc = metadata.encoding.encodingLimits.slice.maximum+1;

        params.ImageSize = [metadata.encoding.reconSpace.matrixSize.x, metadata.encoding.reconSpace.matrixSize.y, metadata.encoding.reconSpace.matrixSize.z];

        params.FOVx_mm = [metadata.encoding.reconSpace.fieldOfView_mm.x, metadata.encoding.reconSpace.fieldOfView_mm.y, metadata.encoding.reconSpace.fieldOfView_mm.z];

        params.PixelSpacing_mm = params.FOVx_mm./params.ImageSize;

        params.NCha = metadata.acquisitionSystemInformation.receiverChannels;

        params.NEco = length(metadata.sequenceParameters.TE);

        params.TE_ms = metadata.sequenceParameters.TE;

        fov_offset_mm = group{1}.head.position(1:2); % check

        params.ImageOrigin_xy(1) = -params.PixelSpacing_mm(1).*(params.ImageSize(1)/2) + fov_offset_mm(1); % mid-pixel - offset
        params.ImageOrigin_xy(2) = -params.PixelSpacing_mm(2).*(params.ImageSize(2)/2) + fov_offset_mm(2); % mid-pixel - offset
      
        % order interleaved slices
        slicePositions = group{1}.head.position(3);

        [~, sliceIndicesSorted] = sort(slicePositions, 'ascend'); % to do: fix for multi

        params.SliceIndicesSorted = sliceIndicesSorted;

        % Format data into a single [RO PE cha] array
        ksp = double(cell2mat(permute(cellfun(@(x) x.data, group, 'UniformOutput', false), [1 3 2])));
        
        % Format navigator data
        navdata = double(cell2mat(permute(cellfun(@(x) x.data, navgroup, 'UniformOutput', false), [1 3 2])));
    
        nColNav = size(navdata,1);
        nNav = length(navgroup);
        sSMP = 130;
        w=10;
        
        navData = reshape(navdata, nColNav, params.NCha, params.NSlc, []);

        navData = navData(:,:,params.SliceIndicesSorted, :);

        fidnav = squeeze(mean(navData(sSMP-w/2:sSMP+w/2,:,:,:))); % cha, slc, lin

        if params.NSlc==1
            fidnav = reshape(fidnav, params.NCha, 1, []);
        end
        
        %{
        figure; plot(abs(squeeze(fidnav(:,1,:))).');

        % Check correlation between simulated and measured FIDnav (should lie along
        % identity line)
        iSli = 1;
        figure; subplot(121); compareFIDs(fidnav(:,iSli,30), calibrationMatrixComplex(:,1,iSli), 'abs');
        subplot(122); compareFIDs(fidnav(:,iSli,30), calibrationMatrixComplex(:,1,iSli), 'angle');
        %}
        
        % To speed up sample every 4th point
        fidnav(:,:,1:25) = repelem(fidnav(:,:,26),1,1,25); % remove points in approach to steady state
        fidnav_smp = fidnav(:,:,2:4:end);
        nT = size(fidnav_smp,3);

        order = 2; nCoeff = 5; 
        
        calibrationMatrix = cat(1, real(calibrationMatrixComplex), imag(calibrationMatrixComplex));

        % Calculate shim parameters using linear model
        fidnav_shim_params = zeros(nT, nCoeff, params.NSlc);

        for iSli=1:params.NSlc

            yy = squeeze(fidnav_smp(:,iSli,:));

            Y = cat(1, real(yy), imag(yy));

            X = lscov(calibrationMatrix(:,:,iSli),Y);

            fidnav_shim_params(:,:,iSli) = X(2:end,:).'; % units uT and uT/m

        end

    % Augmented SENSE recon

        % GenerateMask
        My = params.ImageSize(1); Mx = params.ImageSize(2);  inc = My/nT;

        % Mask k-space where lines were acquired
        shim_mask = zeros(Mx,My,nT); 
        n=1;
        for iT = 1:nT
            shim_mask(:,n:n+inc-1,iT)=1;
            n=n+inc;
        end

        nVC = 16; % set number of virtual channels for recon

        % reshape into slices
        kdata = reshape(ksp, size(ksp,1), params.NCha, params.NEco, params.NSlc, []);

        kdata = permute(kdata, [1 5 4 3 2]);

        kdata = permute(kdata, [permutedims, 4, 5]); 

        if swapdims(1)
            kdata = circshift(flip(kdata,1),1,1);
        end

        if swapdims(2)
            kdata = circshift(flip(kdata,2),1,2);
        end

       
        kdata = kdata(:,:,sliceIndicesSorted,:,:);
        
        imdata = zeros(params.ImageSize(1), size(kdata,2), params.NSlc, params.NEco, params.NCha);

        % Get field maps 
        % Apply same transform to shim mask
        shim_mask = permute(shim_mask,permutedims);

        if swapdims(1)
            shim_mask = circshift(flip(shim_mask,1),1,1);
        end

        if swapdims(2)
            shim_mask = circshift(flip(shim_mask,2),1,2);
        end

        gamma = 42.58; % Hz/uT

        % Get b0maps from shim params
        fprintf('Getting B0 maps from shim parameters...\n');
        b0_map_rad = zeros(params.ImageSize(1), params.ImageSize(2), nT, params.NSlc); 

        for iS = 1:params_cal.NSlc
            for iT = 1:nT
                b0_map_uT = createSphericalHarmonicsFieldMap2d(fidnav_shim_params(iT,:,iS), params.ImageSize(1:2), params.PixelSpacing_mm(1:2), params.ImageOrigin_xy);
                b0_map_rad(:,:,iT,iS) = 2*pi*gamma*b0_map_uT;
            end
        end

        TE_s = params.TE_ms*1e-3;

        imrecon_iter = zeros(params.ImageSize(1), size(kdata,2), params.NSlc, params.NEco, nVC); % params.NEco
        imrecon_ifft = zeros(size(imrecon_iter));

        for iC = 1:nVC

            for iE = 1:params.NEco

                for iS = 1:params.NSlc

                    %imdata(:,:,iZ,iC) = fftshift(ifft2(ifftshift(kdata(1:2:end,:,iZ,iC)))); 
                    b0_aug = zeros(params.ImageSize(1),params.ImageSize(2),1,nT);

                    for iT = 1:nT
                        b0_aug(:,:,:,iT) = ones(params.ImageSize(1),params.ImageSize(2),1).*exp(1i*b0_map_rad(:,:,iT,iS)*TE_s(iE));
                    end

                    if iC == 1

                        thisDataAllC = permute(squeeze(kdata(1:2:end,:,iS,iE,:)), [1 3 2]);              

                            if iE == 1

                                [virtualData, compressionMatrix] = compressCoilData(thisDataAllC,nVC,2);

                            else % apply same compression

                                virtualData = applyCompressionCoilData(thisDataAllC, compressionMatrix, 2);

                            end

                    end % first channel only

                    %thisData2d = squeeze(kdata(1:2:end,:,iS,iE,iC));
                    thisData2d = squeeze(virtualData(:,iC,:,:,:));

                    t_start_recon = tic;
                        imrecon_iter(:,:,iS,iE,iC) = tew_cartesian_augmented_SENSE(double(thisData2d),ones([params.ImageSize(1), params.ImageSize(2)]),double(b0_aug),shim_mask); 
                    t_end_recon = toc(t_start_recon);

                    fprintf('Reconstructed data from Channel %d/%d, Echo %d/%d, Slice %d/%d in %.2f seconds\n', iC, nVC, iE, params.NEco, iS, params.NSlc, t_end_recon);

                    imrecon_ifft(:,:,iS,iE,iC) =  ismrm_transform_kspace_to_image(double(thisData2d),[1,2]); % same recon...

                end

            end

        end
        
        
      % Save iterative recon
      for iEco = 1:params.NEco
      
          imrecon = adapt_array_2d(squeeze(imrecon_iter(:,:,1,iEco,:)), eye(nVC), 0);

          imrecon = abs(imrecon);

          imrecon = imrecon.*(32767./max(imrecon(:)));

          imrecon = int16(round(imrecon));

          image2d = ismrmrd.Image(imrecon);

        % Copy the relevant AcquisitionHeader fields to ImageHeader
        N = params.ImageSize;

        image2d.head = image2d.head.fromAcqHead(group{1}.head);

        image2d.head.data_type = uint16(2);

        image2d.head.image_series_index = 1;

        image2d.head.image_index = iEco;
        
        % Set ISMRMRD Meta Attributes
        meta = struct;
        meta.DataRole               = 'Image';
        meta.ImageProcessingHistory = 'MATLAB';
        meta.WindowCenter           = uint16(16384);
        meta.WindowWidth            = uint16(32768);
        meta.ImageRowDir            = image2d.head.read_dir; 
        meta.ImageColumnDir         = image2d.head.phase_dir;
               
        meta.Keep_image_geometry = {'bool', true}; % this is needed to keep image consistent with labels otherwise flipped
        
        % set_attribute_string also updates attribute_string_len
        image2d = image2d.set_attribute_string(ismrmrd.Meta.serialize(meta));
        
        images{iEco} = image2d;
        
      end
      
      % Save conventional recon
      for iEco = 1:params.NEco
      
          imrecon = adapt_array_2d(squeeze(imrecon_ifft(:,:,1,iEco,:)), eye(nVC), 0);

          imrecon = abs(imrecon);

          imrecon = imrecon.*(32767./max(imrecon(:)));

          imrecon = int16(round(imrecon));

          image2d = ismrmrd.Image(imrecon);

        % Copy the relevant AcquisitionHeader fields to ImageHeader
        N = params.ImageSize;

        image2d.head = image2d.head.fromAcqHead(group{1}.head);

        image2d.head.data_type = uint16(2);

        image2d.head.image_series_index = 2;

        image2d.head.image_index = iEco;
        
        % Set ISMRMRD Meta Attributes
        meta = struct;
        meta.DataRole               = 'Image';
        meta.ImageProcessingHistory = 'MATLAB';
        meta.WindowCenter           = uint16(16384);
        meta.WindowWidth            = uint16(32768);
        meta.ImageRowDir            = image2d.head.read_dir; 
        meta.ImageColumnDir         = image2d.head.phase_dir;
               
        meta.Keep_image_geometry = {'bool', true}; % this is needed to keep image consistent with labels otherwise flipped
        
        % set_attribute_string also updates attribute_string_len
        image2d = image2d.set_attribute_string(ismrmrd.Meta.serialize(meta));
        
        images{iEco+params.NEco} = image2d;
        
      end
             
    end
    
  end
end
