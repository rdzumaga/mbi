function addButton() {
    var btn = document.createElement("BUTTON");
    var t = document.createTextNode("WYKONAJ");
    btn.appendChild(t);
    document.getElementById("buttonPlace").appendChild(btn);
    btn.className = "btn btn-primary btn-lg";
    btn.onclick = function () {
         document.getElementById("seqA").value = document.getElementById("seqA").value.toUpperCase();
         document.getElementById("seqB").value = document.getElementById("seqB").value.toUpperCase();
        if (sequenceValidation()) {
            document.getElementById("myForm").submit();
        }
    };
}

function sequenceValidation() {
    var seqA = document.getElementById("seqA");
    var seqB = document.getElementById("seqB");
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
    if (isNaN(penalty.value) || parseInt(Number(penalty.value)) != penalty.value || isNaN(parseInt(penalty.value, 10))) {
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