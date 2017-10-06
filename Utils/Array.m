classdef Array

    % Class to define an array object. Contains the properties of the
    % electrode array: electrode positions, number of electrodes, electrode
    % spacing, disconnected channels, and sample rate.
    % Lauren Grosberg 6/2016
    properties (SetAccess = private) % Only members of the Array class can set the properties
        xc                   % x coordinates
        yc                   % y coordinates
        numElecs             % number of electrodes
        elecSpacing          % electrode spacing (in microns)
        nullChannels         % disconnected channels
        sampleRate           % sample rate of the array
        difPos               %Distances between electrodes (microns)
        center               %center electrode
        maxR
    end
    
    methods
        function obj = Array(arrayID)
            switch arrayID
                case {1502,1504}
                    [obj.xc,obj.yc] = getElectrodeCoords519();                 
                    obj.numElecs = 519;
                    obj.elecSpacing = 30; % microns
                    obj.nullChannels = [1 130 259 260 389 390 519];
                    obj.sampleRate = 20000;
                    positions=getPositions(obj);
                    x=repmat(nansum(positions'.^2)',1,obj.numElecs);
                    xy=positions*positions';
                    obj.difPos=sqrt(x-2*xy+x');
                    obj.center=448;
                    obj.maxR=nanmax(sqrt(nansum(positions'.^2)));

                case {501,504}

                    [obj.xc,obj.yc] = getElectrodeCoords512();                 
                    obj.numElecs = 512;
                    obj.elecSpacing = 60; % microns
                    obj.nullChannels = [];
                    obj.sampleRate = 20000;
                    positions=getPositions(obj);
                    
                    x=repmat(nansum(positions'.^2)',1,obj.numElecs);
                    xy=positions*positions';
                    obj.difPos=sqrt(x-2*xy+x');
                    obj.center=316;
                    obj.maxR=nanmax(sqrt(nansum(positions'.^2)));
                case 1
                    [obj.xc,obj.yc] = getElectrodeCoords61();                 
                    obj.numElecs = 64;
                    obj.elecSpacing = 60; % microns
                    obj.nullChannels = [9 25 57];
                    obj.sampleRate = 20000;
                    positions=getPositions(obj);
                  
                    x=repmat(nansum(positions'.^2)',1,obj.numElecs);
                    xy=positions*positions';
                  
                    center=41;
                    obj.maxR=nanmax(sqrt(nansum(positions'.^2)));
                    obj.difPos=sqrt(x-2*xy+x');
                    obj.center=41;
            end
        end
        
        function numElecs = getNumElecs(obj)
            % Method to return the number of electrodes in the Array
            numElecs = obj.numElecs; 
        end
        
        function elecSpacing = getSpacing(obj)
            % Method to return the electrode spacing in the Array
            elecSpacing = obj.elecSpacing; 
        end
        
        function positions = getPositions(obj)
            % Method to return the electrode positions
            % size(positions) = numElecs x 2
            positions = [obj.xc; obj.yc]'; 
        end
        
        function activeElectrodes = getElectrodes(obj)
            % Method to return the active electrodes (excludes disconnected
            % channels)
            activeElectrodes = setdiff(1:obj.numElecs,obj.nullChannels); 
        end
        
        function neighborElectrodes = getNeighbors(obj,centerElectrode,orderOfNeighborhood)
            % Method that returns the n=orderOfNeiborhood neighboring 
            % electrodes of a center electrode. Output includes
            % the center electrode and orders the neighbors in numerical
            % order.
            maxDist = obj.elecSpacing*orderOfNeighborhood;
            distances = zeros(1, length(obj.xc));
            
            for ii = 1:length(obj.xc)
                distances(ii) = norm([obj.xc(centerElectrode) - obj.xc(ii), obj.yc(centerElectrode) - obj.yc(ii)]);
            end
            
            % Find the center electrode and all nearest neighbors using
            % radial distance from the center electrode. The value 1.15 is
            % because the electrode spacing is not exact on our specific
            % hardware.
            neighborElectrodes = find(squeeze(distances) < maxDist*1.15); 
           
            % Sort such that the center electrode is first and sorts the
            % rest are in numerical order
            neighborElectrodes(neighborElectrodes == centerElectrode) = [];
            neighborElectrodes = [centerElectrode sort(neighborElectrodes)]; 
        end
    end
end


