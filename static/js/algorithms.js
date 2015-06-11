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
    var errorA = document.getElementById('errorA');
    var errorB = document.getElementById('errorB');
    var div = document.createElement('div');
    if (seq == 'A' && errorA) {
        errorA.innerHTML= "";
    }
    else if(errorB) {
        errorB.innerHTML= "";
    }
}
