classdef fidnav_b0_calib_pa < handle
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
            
           
            % When this criteria is met, run process_raw() on the accumulated
            % data, which returns images that are sent back to the client.
           % if item.head.flagIsSet(item.head.FLAGS.ACQ_LAST_IN_SLICE)
           if item.head.flagIsSet(item.head.FLAGS.ACQ_LAST_IN_MEASUREMENT)
              logging.info("Processing a group of %d k-space data", length(acqGroup))
              images = obj.process_raw(acqGroup, config, metadata, logging);
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
         image = obj.process_raw(acqGroup, config, metadata, logging);
         logging.debug("Sending image to client")
         connection.send_image(image);
         acqGroup = cell(1,0);
       end
      
      connection.send_close();
      return
    end
   
    
    function images = process_raw(obj, group, config, metadata, logging)
      
        % Reconstruct low-resolution PA images
        permutedims = [1 2 3];
        swapdims = [1 0]; % pa
      
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
      
     
        % Format data into a single [RO PE cha] array
        ksp = double(cell2mat(permute(cellfun(@(x) x.data, group, 'UniformOutput', false), [1 3 2])));
        % ksp = permute(ksp, [1 3 2]);
      
        % reshape into slices
        
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

        % order interleaved slices
        slicePositions = group{1}.head.position(3);

        [~, sliceIndicesSorted] = sort(slicePositions, 'ascend'); % to do: fix for multi

        params.SliceIndicesSorted = sliceIndicesSorted;

        kdata = kdata(:,:,sliceIndicesSorted,:,:);
        
        imdata = zeros(params.ImageSize(1), size(kdata,2), params.NSlc, params.NEco, params.NCha);

        for iC = 1:params.NCha

            for iE = 1:params.NEco

                for iZ = 1:params.NSlc

                    imdata(:,:,iZ,iE,iC) = fftshift(ifft2(ifftshift(kdata(1:2:end,:,iZ,iE,iC)))); 

                end

            end
        end

       % imdata = permute(imdata,[2 1 3 4 5]);
       imdata_pa = imdata;
       
       load('/tmp/cal_ap.mat'); % imdata_ap
       
       imdata_pa = permute(imdata_pa,[2 1 3 4 5]);
       
       imdata_ap = permute(imdata, [2 1 3 4 5]);
       
        nSli = params.NSlc;

        Nc = params.NCha;

        iEco = 2;

        
        % Phase correct and average AP/PA calibration images
        phi = angle( sum( (imdata_ap(:,:,:,iEco,:).*conj(imdata_pa(:,:,:,iEco,:))) ,4)/Nc )/2;
        imdata_ap_corr = imdata_ap(:,:,:,iEco,:).*exp(-1i*phi);
        imdata_pa_corr = imdata_pa(:,:,:,iEco,:).*exp(1i*phi);
        coil_im = (imdata_ap_corr + imdata_pa_corr)/2;

        % Create SH basis functions
        % order = 1; nCoeff = 3; % for now just do first order - f0, X, Y
        order = 2; nCoeff = 5; 

        N = params.ImageSize(1:2);

        B0map = zeros([N nCoeff+1]);
        coeffs = cat(1, zeros(1,nCoeff), eye(nCoeff)); % f0 uT, X uT/m, Y uT/m, h.o.t uT/m2

        for ii = 1:nCoeff+1
            B0map(:,:,ii) = createSphericalHarmonicsFieldMap2d(coeffs(ii,:), params.ImageSize(1:2), params.PixelSpacing_mm(1:2), params.ImageOrigin_xy);
        end

        % B0 map is size 64 x 64 x nCoeffs+1, B is vector size Np x nCoeffs+1
        B = reshape(B0map, [], nCoeff+1); % m_sh nPoints*4

        gamma = 42.58*1e6; % Hz/T

        W = 1e-6*B*(2*pi*gamma); % scale (radians)

        t_s = params.TE_ms(iEco)*1e-3; % echo time in seconds

        % Generate calibration matrix for each slice
        for iS = 1:nSli

            Q = reshape(coil_im(:,:,iS,1,:), [], Nc).'; % 64xNp complex vector

            QB{iS} = Q*exp(1i*W*t_s); % 64x4 complex vector

            QBri{iS} = cat(1, real(QB{iS}), imag(QB{iS})); % 128 x 4 concatenate real imaginary

        end

        calibrationMatrix = cat(3, QBri{:});
        calibrationMatrixComplex = complex(calibrationMatrix(1:Nc,:,:), calibrationMatrix(Nc+1:end,:,:));

        % Save calibration matrix to tmp directory
        save('/tmp/calibrationMatrix.mat', 'calibrationMatrixComplex', 'params');
        
 %       imdata = permute(imdata, [2 1 3 4]); % because Matlab co-ordinates are Y-X-Z
      iEco = 2;
      
      imrecon = adapt_array_2d(squeeze(imdata_pa(:,:,1,iEco,:)), eye(params.NCha), 0);
      
      imrecon = permute(abs(imrecon),[2 1]);
      
      imrecon = imrecon.*(32767./max(imrecon(:)));
            
      imrecon = int16(round(imrecon));

      image2d = ismrmrd.Image(imrecon);
        
    % Copy the relevant AcquisitionHeader fields to ImageHeader
    N = params.ImageSize;
    
    image2d.head = image2d.head.fromAcqHead(group{1}.head);

    image2d.head.data_type = uint16(2);

    % image2d.head.image_series_index = 2;

    % image2d.head.image_index = 1;
        
        % field_of_view is mandatory
   %{
     image2d.head.field_of_view  = single([metadata.encoding(1).reconSpace.fieldOfView_mm.x ...
          metadata.encoding(1).reconSpace.fieldOfView_mm.y ...
          metadata.encoding(1).reconSpace.fieldOfView_mm.z]);
        
        image2d.head.matrix_size(1) = uint16(N(1));
        image2d.head.matrix_size(2) = uint16(N(2));
        image2d.head.matrix_size(3) = uint16(1);
     %}   
        % THIS REPLICATES ICE RECON BEHAVIOR WITH NORM ORIENT OFF
       % image2d.head.read_dir = group{1}.head.read_dir; % removed -1
       % image2d.head.phase_dir = group{1}.head.phase_dir; % removed -1
                           
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
        
        images{1} = image2d;
             
    end
    
  end
end
