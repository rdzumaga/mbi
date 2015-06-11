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
    if (!seqA.value.match("A|C|T|G")) {
        seqA.className = "nonValidate";
        addWarning('A');
        ok = false;
    }
    if (!seqB.value.match("A|C|T|G")) {
        seqB.className = "nonValidate";
        addWarning('B');
        ok = false;
    }
    return ok;
}

function addWarning(seq) {
    var div = document.createElement('div');
    div.className = 'error';
    if (seq == 'A') {
        div.innerHTML = '<div class="error" id="break_penalty__error" style="display: inline-block;">Sekwencja A jest nieprawidłowa! Musi składać się wyłącznie ze znaków A,C,T,G i nie może być pusta</div>';
        document.getElementById('seqADiv').appendChild(div);
    }
    else {
        div.innerHTML = '<div class="error" id="break_penalty__error" style="display: inline-block;">Sekwencja B jest nieprawidłowa! Musi składać się wyłącznie ze znaków A,C,T,G i nie może być pusta</div>';
        document.getElementById('seqBDiv').appendChild(div);
    }
}

function removeWarning(input) {
    document.getElementById('content').removeChild(input.parentNode);
}
