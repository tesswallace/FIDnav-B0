classdef server < handle

    properties
        port      = [];
        %tcpHandle = [];
        server_socket = [];
        output_socket = [];
        log       = [];
    end

    methods
        function obj = server(port,log)
            log.info('Initializing server on port %d', port);
            obj.port = port;
            obj.log = log;
        end

        function serve(obj)
            while true
                try
                    obj.server_socket = ServerSocket(obj.port);
                    obj.log.info('Waiting for client to connect to this host on port : %d', obj.port);
                    [obj.output_socket,remote_socket_address] = SocketAccept(obj.server_socket);
                    obj.log.info('Accepting connection from: %s', remote_socket_address);
                    handle(obj);
                    SocketClose(obj.server_socket);
                    SocketClose(obj.output_socket);
                catch
                    if ~isempty(obj.output_socket)
                        SocketClose(obj.output_socket);
                    end
                    if ~isempty(obj.server_socket)
                        SocketClose(obj.server_socket);
                    end
                    %pause(1);
                end
                pause(1)
            end
        end

        function handle(obj)
            try
                % conn = connection(obj.tcpHandle, obj.log);
                conn = connection(obj.output_socket, obj.log);
                config = next(conn);
                metadata = next(conn);

                try
                    metadata = ismrmrd.xml.deserialize(metadata);
                    if ~isempty(metadata.acquisitionSystemInformation.systemFieldStrength_T)
                        obj.log.info("Data is from a %s %s at %1.1fT", metadata.acquisitionSystemInformation.systemVendor, metadata.acquisitionSystemInformation.systemModel, metadata.acquisitionSystemInformation.systemFieldStrength_T)
                    end
                catch
                    obj.log.info("Metadata is not a valid MRD XML structure.  Passing on metadata as text")
                end

                % Decide what program to use based on config
                % As a shortcut, we accept the file name as text too.
                if strcmpi(config, "simplefft")
                    obj.log.info("Starting simplefft processing based on config")
                    recon = simplefft;
                elseif strcmpi(config, "invertcontrast")
                    obj.log.info("Starting invertcontrast processing based on config")
                    recon = invertcontrast;
                elseif strcmpi(config, "mapvbvd")
                    obj.log.info("Starting mapvbvd processing based on config")
                    recon = fire_mapVBVD;
                    %{
                elseif strcmpi(config, "reconradialtew")
                    obj.log.info("Starting reconradialtew processing based on config")
                    recon = reconradialtew;
                elseif strcmpi(config, "reconradialtfnufft")
                    obj.log.info("Starting reconradialtfnufft processing based on config")
                    recon = reconradialtfnufft;
                elseif strcmpi(config, "reconradialtfnufftv2")
                    obj.log.info("Starting reconradialtfnufftv2 processing based on config")
                    recon = reconradialtfnufftv2;                    
                elseif strcmpi(config, "reconradialtfnufftstream")
                    obj.log.info("Starting reconradialtfnufftstream processing based on config")
                    recon = reconradialtfnufftstream;
                elseif strcmpi(config, "reconradialtfnufftstreamv2")
                    obj.log.info("Starting reconradialtfnufftstreamv2 processing based on config")
                    recon = reconradialtfnufftstreamv2;
                elseif strcmpi(config, "reconradialtfnufftstreamv2_moco")
                    obj.log.info("Starting reconradialtfnufftstreamv2_moco processing based on config")
                    recon = reconradialtfnufftstreamv2_moco;   
                elseif strcmpi(config, "fidnav_kidney_calibration")
                    obj.log.info("Starting fidnav_kidney_calibration processing based on config")
                    recon = fidnav_kidney_calibration;
                %}
                elseif strcmpi(config, "fidnav_b0_calib_ap")
                    obj.log.info("Starting fidnav_b0_calib_ap processing based on config")
                    recon = fidnav_b0_calib_ap;   
                elseif strcmpi(config, "fidnav_b0_calib_pa")
                    obj.log.info("Starting fidnav_b0_calib_pa processing based on config")
                    recon = fidnav_b0_calib_pa; 
                elseif strcmpi(config, "fidnav_b0_retro_recon")
                    obj.log.info("Starting fidnav_b0_retro_recon processing based on config")
                    recon = fidnav_b0_retro_recon;     
                else
                    if exist(config, 'class')
                        obj.log.info("Starting %s processing based on config", config)
                        eval(['recon = ' config ';'])
                    else
                        obj.log.info("Unknown config '%s'.  Falling back to 'invertcontrast'", config)
                        recon = invertcontrast;
                    end
                end
                recon.process(conn, config, metadata, obj.log);

            catch ME
                obj.log.error('[%s:%d] %s', ME.stack(2).name, ME.stack(2).line, ME.message);
                rethrow(ME);
            end
        end

        function delete(obj)
           if ~isempty(obj.server_socket)
                SocketClose(obj.server_socket);
            end
        end

    end

end
