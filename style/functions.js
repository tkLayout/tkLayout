
function pippo() {
    alert("pippo");
}

// DO NOT USE THIS
function getElementbyClassName(container, className) {
    var result=new Array();
    var iFound=0;
    for (var iElement=0; iElement<container.length; iElement++) {
	if (container[iElement].className==className)
	    result[iFound++]=container[iElement];
	var subFound = new Array();
	subFound = getElementbyClassName(container[iElement], className);
	for (var iSubElement=0; iSubElement<subFound.length; iSubElement++) {
	    result[iFound++] = subFound[iSubElement];
	}
    }
    return result;
}

function returnElementById( elementId ) {
    var result;
    if (document.getElementById)
        var result = document.getElementById(elementId);
    else if (document.all)
        var result = document.all[elementId];
    else if (document.layers)
        var result = document.layers[elementId];
    return result;
}

function setRowOver(evt) {
    this.className = "rowOver";
}

function setRowA(evt) {
    this.className = "rowA";
}

function setRowB(evt) {
    this.className = "rowB";
}

function stripeTable(myTable) {
    // Get a pointer to the table
    //var myTable = returnElementById(elementId);
    //if (! myTable) { return; }

    
    // Go through the list of tbody
    var myBodies = myTable.getElementsByTagName("tbody");
    for (var iBody = 0; iBody < myBodies.length; iBody++) {
	// Go throught the rows of the body
	var myRows = myBodies[iBody].getElementsByTagName("tr");
	for (var iRow = 0; iRow < myRows.length; iRow++) {
	    if (iRow%2) {
		myRows[iRow].className = "rowA";
		//myRows[iRow].onmouseout = setRowA;
	    } else {
		myRows[iRow].className = "rowB";
		//myRows[iRow].onmouseout = setRowB;
	    }
	    //myRows[iRow].onmouseover = setRowOver;
        }
    }
}

function getNextSibling(startObject) {
    var myObject = startObject.nextSibling;
    while (myObject.nodeType!=1) {
	myObject = myObject.nextSibling;
    }
    return myObject;
}

function toggleDisplay(content) {
    var styleSet = (content.style.display!="");
    var isHidden;
    if (styleSet) isHidden = (content.style.display=="none");
    if (!styleSet) isHidden = (content.className=="hidden");
    content.style.display = (isHidden ? "block" : "none");
    //var chevron = content.parentNode
    //    .firstChild.childNodes[1].childNodes[0];
    //chevron.src = chevron.src
    //    .split(expand ? "expand.gif" : "collapse.gif")
    //    .join(expand ? "collapse.gif" : "expand.gif");
}

function toggleNextDisplay() {
    var nextObject = getNextSibling(this);
    toggleDisplay(nextObject);
}

function getWindowWidth() {
    var x = 0;
    if (self.innerHeight) {
	x = self.innerWidth;
    } else if (document.documentElement && document.documentElement.clientHeight) {
	x = document.documentElement.clientWidth;
    } else if (document.body) {
	x = document.body.clientWidth;
    }
    return x;
}

function getWindowHeight() {
    var y = 0;
    if (self.innerHeight) {
	y = self.innerHeight;
    } else if (document.documentElement && document.documentElement.clientHeight) {
	y = document.documentElement.clientHeight;
    } else if (document.body) {
	y = document.body.clientHeight;
    }
    return y;
}

function hidePopup(evt) {
    var myObject;
    myObject = returnElementById("popupDiv");
    if (myObject) {
	myObject.style.display = "none";
    }
    myObject = returnElementById("curtainDiv");
    if (myObject) {
	myObject.style.display = "none";
    }
}

function setOpacity(myObject, value) {
    myObject.style.opacity = value;
    myObject.style.filter = 'alpha(opacity=' + value + ')';
}


function animateFade(objectId, endTime, startValue, endValue, timeToFade) {
    var currentTime = new Date().getTime();
    var animationFractionRemaining = (endTime - currentTime) / timeToFade;
    if (animationFractionRemaining<0) animationFractionRemaining=0;
    var myObject = returnElementById(objectId);
    var myOpacity = (startValue - endValue) * animationFractionRemaining + endValue;
    if (myObject) {
	setOpacity(myObject, myOpacity);
	if (animationFractionRemaining!=0)
	    setTimeout("animateFade('" + myObject.id + "',"+endTime+", "+startValue+", "+endValue+", "+timeToFade+")", 30);
    }
}

function fadeOpacity(myObject, startValue, endValue, timeToFade) {
    var currentTime = new Date().getTime();
    setOpacity(myObject, startValue);
    animateFade(myObject.id, currentTime+timeToFade, startValue, endValue, timeToFade);
}

function curtainIn() {
    var curtainDiv = returnElementById("curtainDiv");
    if (!curtainDiv) {
        curtainDiv=document.createElement("DIV");
        curtainDiv.id="curtainDiv";
	curtainDiv.className = "curtain";
	curtainDiv.style.display = "none";
	curtainDiv.onclick=hidePopup;
	document.body.appendChild(curtainDiv);
    } else {
	curtainDiv.style.display = "none";
    }
    curtainDiv.style.left="0px";
    curtainDiv.style.top="0px";
    //curtainDiv.style.width="100%";
    //curtainDiv.style.height="100%";
    curtainDiv.style.height=document.body.scrollHeight;
    curtainDiv.style.width=document.body.scrollWidth;
    curtainDiv.style.position="absolute";
    curtainDiv.style.background='black';
    setOpacity(curtainDiv, 0);
    //setOpacity(curtainDiv, 0.6);
    curtainDiv.style.display = "block";
    fadeOpacity(curtainDiv, 0, 0.7, 500);
    return curtainDiv;
}

function addFileLink(baseNode, URLbase, extension, pretext) {
    // Create a link for the plot image
    var link = document.createElement('a');
    link.setAttribute('href', URLbase+"."+extension);
    link.setAttribute('name', 'imgDownloadLink'+URLbase+"."+extension);
    baseNode.appendChild(document.createTextNode(pretext));
    baseNode.appendChild(link);
    // Put the text in the link
    link.appendChild(document.createTextNode("("+extension+")"));
}

function popupDiv(URLbase, myWidth, myHeight, title, extFile) {
    var curtainDiv = curtainIn();
    var popupDiv = returnElementById("popupDiv");
    if (!popupDiv) {
        popupDiv = document.createElement("DIV");
        popupDiv.id = "popupDiv";
	popupDiv.style.display = "none";
	//popupDiv.onclick = hidePopup;
	document.body.appendChild(popupDiv);
    } else {
	popupDiv.style.display = "none";
    }

    var popupImg = returnElementById("popupImg");
    if (!popupImg) {
	popupImg = document.createElement("IMG");
	popupImg.id = "popupImg";
	//popupImg.onclick = hidePopup;
	popupDiv.appendChild(popupImg);
    }

    var titleDiv = returnElementById("titleDiv");
    if (titleDiv) {
	//removeChildren(titleDiv);
	titleDiv.parentNode.removeChild(titleDiv);
    }

    // Create a DIV for the plot title and links
    titleDiv = document.createElement("DIV");
    titleDiv.id="titleDiv";
    titleDiv.style.clear="left";
    popupDiv.appendChild(titleDiv);	

    // Put the title in the title div, if exists
    if (title) {
        var myText = document.createTextNode(title+" - ");
        titleDiv.appendChild(myText);
    }

    // Add links to filetypes (the png is there by default)
    addFileLink(titleDiv, URLbase, "png", "");
    if (extFile) { 
      var x;
      var myTypes = new Array();
      myTypes = extFile.split("|");
      for (x in myTypes) {
        addFileLink(titleDiv, URLbase, myTypes[x], " - ");
      }
    }

    var myLeft = Math.floor((getWindowWidth()-myWidth)/2);
    var myTop = Math.floor((getWindowHeight()-myHeight)/2)+document.body.scrollTop;
    popupDiv.style.left=myLeft+"px";
    popupDiv.style.top=myTop+"px";
    //popupDiv.style.width=myWidth+"px";
    //popupDiv.style.height=myHeight+"px";
    popupDiv.style.padding="10px";

    popupImg.style.width=popupDiv.style.width;
    popupImg.style.height=popupDiv.style.height;

    popupDiv.style.position="absolute";
    //popupDiv.style.background='url('+URL+')';
    popupDiv.style.background='white';
    popupImg.src=URLbase+".png";
    popupImg.border=0;

    popupDiv.style.display = "block";
}

function preparePageGoodies() {
    // Make the tables striiped
    var allTables = document.getElementsByTagName('table');
    for (var iTable = 0; iTable < allTables.length; iTable++) {
	stripeTable(allTables[iTable]);
    }

    // Expanding h2 titles
    var allTitles = document.getElementsByTagName('h2');
    for (var iTitle = 0; iTitle < allTitles.length; iTitle++) {
	if  (allTitles[iTitle].className=="hidingTitle") {
	    allTitles[iTitle].onclick = toggleNextDisplay;
	}
    }

}
