#!/bin/bash

if [ $# != 1 ]; then
  echo Syntax: $0 materialName
  exit -1
fi

materialName=$1
fileName=propagate_$materialName

if [ -f $fileName ]; then
  echo $fileName already exists
  exit -1
fi

cat > $fileName << @EOF
// 1:1 conversion of $materialName

Conversion {
  Input {
    Element {
      elementName $materialName 
      quantity 1
      unit g/m
    }
  }
  Output {
    Element {
      elementName $materialName
      quantity 1
      unit g/m
      service true
    }
  }
}
@EOF
