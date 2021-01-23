jQuery('button').click( function(e) {
    jQuery('.collapse').collapse('hide');
});

const formSubmit = document.querySelector("#submission-form");
const pdbFile = document.querySelector('#pdb');
const Smiles = document.querySelector('#smiles');
const submitButton = document.querySelector('#submit-button');
const msg = document.querySelector('.msg');
const Spinner = document.querySelector('#spinner');

function errorMessage() {
    msg.innerHTML = "Provide either a valid PDB or SMILES string of your molecule.";
    msg.style.color = "red";
    msg.style.paddingTop = "10px";
    setTimeout(function() {
        msg.innerHTML=null;
        msg.style.paddingTop= "0px";
    }, 3000);

};

function showSpinner() {

    submitButton.style.visibility = "hidden";
    Spinner.style.visibility = "visible";

};

formSubmit.addEventListener('submit', function(e) {

    if (Smiles.value) {
        showSpinner();
    } else if (pdbFile.value.split('.').pop() == "pdb") {
        showSpinner();
    } else {
        e.preventDefault();
        errorMessage();
    }

});

