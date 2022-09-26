function [box]= Boxing (vidFrame, yIndex,xIndex,iBox)



    box = vidFrame(yIndex-iBox:yIndex+iBox,xIndex-iBox:xIndex+iBox);
