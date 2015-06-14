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
    var seqLength = document.getElementById("seqLength");
    var penalty = document.getElementById("break_penalty");
    var ok = true;
    if (seqA.value=="" || !(/^[ACTG]*$/.test(seqA.value))) {
        addWarning('A');
        ok = false;
    }
    else {
        removeWarning('A');
    }
    if (seqB.value=="" || !(/^[ACTG]*$/.test(seqB.value))) {
        addWarning('B');
        ok = false;
    }
    else {
        removeWarning('B');
    }
    if (isNaN(seqLengthseqLength.value) || parseInt(Number(seqLength.value)) != seqLength.value || isNaN(parseInt(seqLength.value, 10)) || seqLength.value <1 || seqLength.value >6) {
        addWarning('C');
        ok = false;
    }
    else {
        removeWarning('C');
    }
    return ok;
}

function addWarning(seq) {
    var errorA = document.getElementById('errorA');
    var errorB = document.getElementById('errorB');
    var errorC = document.getElementById('errorC');
    
    if (seq == 'A'){
        if (!errorA) {
            var div = document.createElement('div');
            div.innerHTML = '<div class="error" id="errorA" style="display: inline-block;">Sekwencja A jest nieprawidłowa! Musi składać się wyłącznie ze znaków A,C,T,G i nie może być pusta!</div>';
            document.getElementById('seqADiv').appendChild(div);
        }
        else {
            errorA.innerHTML = "Sekwencja A jest nieprawidłowa! Musi składać się wyłącznie ze znaków A,C,T,G i nie może być pusta!";
        }
    }
    else if (seq == 'B'){
        if (!errorB) {
            var div = document.createElement('div');
            div.innerHTML = '<div class="error" id="errorB" style="display: inline-block;">Sekwencja B jest nieprawidłowa! Musi składać się wyłącznie ze znaków A,C,T,G i nie może być pusta!</div>';
            document.getElementById('seqBDiv').appendChild(div);
        }
        else {
            errorB.innerHTML = "Sekwencja B jest nieprawidłowa! Musi składać się wyłącznie ze znaków A,C,T,G i nie może być pusta!";
        }
    }
     else {
        if (!errorC) {
            var div = document.createElement('div');
            div.innerHTML = '<div class="error" id="errorC" style="display: inline-block;">Wartość kary musi być ujemną liczbą całkowitą!</div>';
            document.getElementById('seqCDiv').appendChild(div);
        }
        else {
            errorC.innerHTML = "Wartość kary musi być ujemną liczbą całkowitą!";
        }
     }
}

function removeWarning(seq) {
    var errorA = document.getElementById('errorA');
    var errorB = document.getElementById('errorB');
    var errorC = document.getElementById('errorC');

    if (seq == 'A' && errorA) {
        errorA.innerHTML= "";
    }
    else if(seq == 'B' && errorB) {
        errorB.innerHTML= "";
    }
    else if (errorC){
        errorC.innerHTML="";
    }
}
