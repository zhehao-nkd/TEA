classdef (Abstract) EphysAnalysis < handle
    properties
        
        names% include  birdname, neuronname, experimentname, channelname, unitname, stimuliname
        dataset %The core data, include stimuli, corresponing response, firing rates,and all other calculations
        
        
        
    end
    
    
     methods
         
         add(obj)
         
         Three(obj) % draw thespectrogram, raster and spike density function
         collectImages(obj) % collect images(raster plots or others) to struct and pass to the next level
         judgeResponse(obj) % judge whether the tested neuron respond to the stimuli or not
         
         
     end
  
end