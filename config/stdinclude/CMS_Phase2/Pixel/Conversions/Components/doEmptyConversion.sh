#!/bin/bash

materialName="$1"

cat << @EOF
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
