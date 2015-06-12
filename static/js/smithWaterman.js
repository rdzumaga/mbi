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
         document.getElementById("seqA").value = document.getElementById("seqA").value.toUpperCase();
         document.getElementById("seqB").value = document.getElementById("seqB").value.toUpperCase();
        if (sequenceValidation()) {
            checkStepMode();
            document.getElementById("myForm").submit();
        }
    };
}

function sequenceValidation() {
    var seqA = document.getElementById("seqA");
    var seqB = document.getElementById("seqB");
    var penalty = document.getElementById("break_penalty");
    var match = document.getElementById("match");
    var mismatch = document.getElementById("mismatch");
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
    if (isNaN(penalty.value) || parseInt(Number(penalty.value)) != penalty.value || isNaN(parseInt(penalty.value, 10)) || penalty.value >=0) {
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
    if (isNaN(mismatch.value) || parseInt(Number(mismatch.value)) != mismatch.value || isNaN(parseInt(mismatch.value, 10)) || mismatch.value >0) {
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
     else if (seq == 'C') {
        if (!errorC) {
            var div = document.createElement('div');
            div.innerHTML = '<div class="error" id="errorC" style="display: inline-block;">Wartość kary musi być ujemną liczbą całkowitą!</div>';
            document.getElementById('seqCDiv').appendChild(div);
        }
        else {
            errorC.innerHTML = "Wartość kary musi być ujemną liczbą całkowitą!";
        }
     }
    else if (seq == 'D') {
        if (!errorD) {
            var div = document.createElement('div');
            div.innerHTML = '<div class="error" id="errorD" style="display: inline-block;">Wartość zgodności musi być dodatnią liczbą całkowitą!</div>';
            document.getElementById('seqDDiv').appendChild(div);
        }
        else {
            errorD.innerHTML = "Wartość zgodności musi być dodatnią liczbą całkowitą!";
        }
     }
    else {
        if (!errorE) {
            var div = document.createElement('div');
            div.innerHTML = '<div class="error" id="errorE" style="display: inline-block;">Wartość kary musi być liczbą całkowitą <=0 !</div>';
            document.getElementById('seqEDiv').appendChild(div);
        }
        else {
            errorE.innerHTML = "Wartość kary musi być liczbą całkowitą <=0!";
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
    else if (seq == 'C' && errorC){
        errorC.innerHTML="";
    }
    else if (seq == 'D' && errorD){
        errorD.innerHTML="";
    }
    else if (seq == 'E' && errorE){
        errorE.innerHTML="";
    }
}

function checkStepMode() {
    if(document.getElementById("checkBox").checked)
    {
        if (document.getElementById("iters").value == -1){
            document.getElementById("iters").value=1;
        }
        else{
            var iters = parseInt(document.getElementById("iters").value) + 1;
            document.getElementById("iters").value = iters;
        }
    }
    else document.getElementById("iters").value=-1;
}


function addValueToMatrix(){
    var iters = document.getElementById("iters").value;
    for(i=1; i<=iters; i++){
        var id = "values" + i;
        var val = document.getElementById(id).value;
        var arrayVal = getArray(val);
        var cellId = "cell" + i;
         document.getElementById(cellId).innerHTML = addArrow(arrayVal[0]) + checkIfStrong(arrayVal, 3) + "<br/>" + checkIfStrong(arrayVal, 4) +  "<br/>" + checkIfStrong(arrayVal, 5);
    }
}

function getArray(val){
    var str = val.substring(1, val.length-1);
    var array = str.split(", ");
    return array;
}

function checkIfStrong(array, nr)
{
    if (array[0] < 3 && (parseInt(array[0])+3) == nr ){
        return "<strong>"+array[nr]+"</strong>";
    }
    else return array[nr];
}

function addArrow (nr){
    switch (nr) {
        case "0":
            return '<i class="glyphicon glyphicon-arrow-right down-right"></i>';
            break;
        case "1":
            return '<i class="glyphicon glyphicon-arrow-down down"></i>';
            break;
        case "2":
            return '<i class="glyphicon glyphicon-arrow-right right"></i>';
            break;
        case "3":
            return '';
            break;
    }
}
