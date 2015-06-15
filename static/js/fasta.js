function addButton() {
    var btn = document.createElement("BUTTON");
    var t = document.createTextNode("WYKONAJ");
    btn.appendChild(t);
    document.getElementById("buttonPlace").appendChild(btn);
    btn.className = "btn btn-primary btn-lg";
    showMaxResult();
    btn.onclick = function () {
         document.getElementById("seq").value = document.getElementById("seq").value.toUpperCase();
        if (sequenceValidation()) {
            document.getElementById("myForm").submit();
        }        
    };
}


function sequenceValidation() {
    
    var seq = document.getElementById("seq");
    var seqLength = document.getElementById("length");
    var thrC = document.getElementById("threshold_c");
    var match = document.getElementById("match");    
    var penalty = document.getElementById("break_penalty");
    
    var ok = true;
    
    if (seq.value=="" || !(/^[ACTG]*$/.test(seq.value))) {
        addWarning('A');
        ok = false;
    }
    else {
        removeWarning('A');
    }
    if (isNaN(seqLength.value) || parseInt(Number(seqLength.value)) != seqLength.value || isNaN(parseInt(seqLength.value, 10)) || seqLength.value <1 || seqLength.value >6) {
        addWarning('B');
        ok = false;
    }
    else {
        removeWarning('B');
    }    
    if (isNaN(thrC.value) || parseInt(Number(thrC.value)) != thrC.value || isNaN(parseInt(thrC.value, 10)) || thrC.value <=0) {
        addWarning('C');
        ok = false;
    }
    else {
        removeWarning('C');
    }    
    if (isNaN(match.value) || parseInt(Number(match.value)) != match.value || isNaN(parseInt(match.value, 10)) || match.value <=0) {
        addWarning('D');
        ok = false;
    }
    else {
        removeWarning('D');
    }
    if (isNaN(penalty.value) || parseInt(Number(penalty.value)) != penalty.value || isNaN(parseInt(penalty.value, 10)) || penalty.value >=0) {
        addWarning('E');
        ok = false;
    }
    else {
        removeWarning('E');
    }
    return ok;
}

function addWarning(seq) {
    
    var errorA = document.getElementById('errorA');
    var errorB = document.getElementById('errorB');
    var errorC = document.getElementById('errorC');
    var errorD = document.getElementById('errorD');
    var errorE = document.getElementById('errorE');
    
    if (seq == 'A'){
        if (!errorA) {
            var div = document.createElement('div');
            div.innerHTML = '<div class="error" id="errorA" style="display: inline-block;">Sekwencja jest nieprawidłowa! Musi składać się wyłącznie ze znaków A,C,T,G i nie może być pusta!</div>';
            document.getElementById('seqDiv').appendChild(div);
        }
        else {
            errorA.innerHTML = "Sekwencja A jest nieprawidłowa! Musi składać się wyłącznie ze znaków A,C,T,G i nie może być pusta!";
        }
    }
    else if (seq == 'B'){
        if (!errorB) {
            var div = document.createElement('div');
            div.innerHTML = '<div class="error" id="errorB" style="display: inline-block;">Długość sekwencji jest nieprawidłowa! Musi być dodatnią liczbą całkowitą z przedziału <1, 6>!</div>';
            document.getElementById('lengthDiv').appendChild(div);
        }
        else {
            errorB.innerHTML = "Długość sekwencji jest nieprawidłowa! Musi być dodatnią liczbą całkowitą z przedziału <1, 6>!";
        }
    }
     else if (seq == 'C'){
        if (!errorC) {
            var div = document.createElement('div');
            div.innerHTML = '<div class="error" id="errorC" style="display: inline-block;">Wartość progu musi być dodatnią liczbą całkowitą!</div>';
            document.getElementById('threshold_cDiv').appendChild(div);
        }
        else {
            errorC.innerHTML = "Wartość progu musi być dodatnią liczbą całkowitą!";
        }
     }
     else if (seq == 'D'){
        if (!errorD) {
            var div = document.createElement('div');
            div.innerHTML = '<div class="error" id="errorD" style="display: inline-block;">Wartość dopasowania musi być dodatnią liczbą całkowitą!</div>';
            document.getElementById('matchDiv').appendChild(div);
        }
        else {
            errorD.innerHTML = "Wartość dopasowania musi być dodatnią liczbą całkowitą!";
        }
     }
    else {
        if (!errorE) {
            var div = document.createElement('div');
            div.innerHTML = '<div class="error" id="errorE" style="display: inline-block;">Wartość kary musi być ujemną liczbą całkowitą!</div>';
            document.getElementById('break_penaltyDiv').appendChild(div);
        }
        else {
            errorE.innerHTML = "Wartość kary musi być ujemną liczbą całkowitą!";
        }
     }
}

function removeWarning(seq) {
    var errorA = document.getElementById('errorA');
    var errorB = document.getElementById('errorB');
    var errorC = document.getElementById('errorC');
    var errorD = document.getElementById('errorD');
    var errorE = document.getElementById('errorE');

    if (seq == 'A' && errorA) {
        errorA.innerHTML= "";
    }
    else if(seq == 'B' && errorB) {
        errorB.innerHTML= "";
    }
    else if(seq == 'C' && errorC) {
        errorC.innerHTML= "";
    }
    else if(seq == 'D' && errorD) {
        errorD.innerHTML= "";
    }
    else if (errorE){
        errorE.innerHTML="";
    }
}

function showMaxResult(){
    var idList = findMaxResult();
    if(idList){
        for	(index = 0; index < idList.length; index++) {
                var id = "targetSpan"+idList[index];
                document.getElementById(id).className = "bestResult";
        }
    }
}

function findMaxResult() {
    var count =  document.getElementById("count").value;
    var tempVal = 0;
    var tempValId = [];
    for(i=0; i<parseInt(count); ++i){
    var id = "result" + i;
     var value = document.getElementById(id);
        if (value){
            var val = value.innerText;
            if (parseInt(val)>=tempVal){

                if (parseInt(val)>tempVal){
                    tempValId.length = 0;
                }
                tempVal = parseInt(val);
                tempValId.push(i);
            }
        }
    }
    return tempValId;
}
