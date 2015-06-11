function addButton() {
    var btn = document.createElement("BUTTON");
    var t = document.createTextNode("WYKONAJ");
    btn.appendChild(t);
    document.getElementById("buttonPlace").appendChild(btn);
    btn.className = "btn btn-primary btn-lg";
    btn.onclick = function () {
        if (sequenceValidation()) {
            document.getElementById("myForm").submit();
        }
    };
}



function sequenceValidation() {
    var seqA = document.getElementById("seqA");
    var seqB = document.getElementById("seqB");
    var ok = true;
    if (!seqA.value.toUpperCase().match("ACTG")) {
        seqA.className = "form-control nonValidate";
        addWarning('A');
        ok = false;
    }
    else {
        removeWarning('A');
    }
    if (!seqB.value.toUpperCase().match("ACTG")) {
        seqB.className = "form-control nonValidate";
        addWarning('B');
        ok = false;
    }
    else {
        removeWarning('B');
    }
    return ok;
}

function addWarning(seq) {
    var errorA = document.getElementById('errorA');
    var errorB = document.getElementById('errorB');
    if (seq == 'A'){
        if (!errorA) {
            var div = document.createElement('div');
            div.innerHTML = '<div class="error" id="errorA" style="display: inline-block;">Sekwencja A jest nieprawidłowa! Musi składać się wyłącznie ze znaków A,C,T,G i nie może być pusta</div>';
            document.getElementById('seqADiv').appendChild(div);
        }
        else {
            errorA.innerHTML = "Sekwencja A jest nieprawidłowa! Musi składać się wyłącznie ze znaków A,C,T,G i nie może być pusta";
        }
    }
    else {
        if (!errorB) {
            var div = document.createElement('div');
            div.innerHTML = '<div class="error" id="errorB" style="display: inline-block;">Sekwencja B jest nieprawidłowa! Musi składać się wyłącznie ze znaków A,C,T,G i nie może być pusta</div>';
            document.getElementById('seqBDiv').appendChild(div);
        }
        else {
            errorB.innerHTML = "Sekwencja B jest nieprawidłowa! Musi składać się wyłącznie ze znaków A,C,T,G i nie może być pusta";
        }
    }
}

function removeWarning(seq) {
    var div = document.createElement('div');
    if (seq == 'A') {
        document.getElementById('errorA').innerHTML= "";
    }
    else {
        document.getElementById('errorB').innerHTML= "";
    }
}
