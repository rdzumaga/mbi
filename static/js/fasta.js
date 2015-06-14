function addButton() {
    var btn = document.createElement("BUTTON");
    var t = document.createTextNode("WYKONAJ");
    btn.appendChild(t);
    document.getElementById("buttonPlace").appendChild(btn);
    btn.className = "btn btn-primary btn-lg";
    
    if (document.getElementById("values1")){
        addValueToMatrix();
    }
    
    btn.onclick = function () {
         document.getElementById("seq").value = document.getElementById("seq").value.toUpperCase();
        if (sequenceValidation()) {
            checkStepMode();
            document.getElementById("myForm").submit();
        }
    };
}


function sequenceValidation() {
    
    var seq = document.getElementById("seq");
    var seqLength = document.getElementById("length");
    var thrC = document.getElementById("threshold_c");
    var thrT = document.getElementById("threshold_t");    
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
    if (isNaN(thrT.value) || parseInt(Number(thrT.value)) != thrT.value || isNaN(parseInt(thrT.value, 10)) || thrT.value <=0) {
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
            div.innerHTML = '<div class="error" id="errorD" style="display: inline-block;">Wartość progu musi być dodatnią liczbą całkowitą!</div>';
            document.getElementById('threshold_tDiv').appendChild(div);
        }
        else {
            errorD.innerHTML = "Wartość progu musi być dodatnią liczbą całkowitą!";
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
