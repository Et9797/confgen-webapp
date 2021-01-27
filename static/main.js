jQuery('button').click( function(e) {
    jQuery('.collapse').collapse('hide');
});

const formSubmit = document.querySelector("#submission-form");
const MolFile = document.querySelector('#mol_file');
const Smiles = document.querySelector('#smiles');
const submitButton = document.querySelector('#submit-button');
const msg = document.querySelector('.msg');
const Spinner = document.querySelector('#spinner');

function errorMessage() {
    msg.innerHTML = "Provide either a valid file format or SMILES string of your molecule.";
    msg.style.color = "red";
    msg.style.paddingTop = "10px";
    setTimeout(function() {
        msg.innerHTML=null;
        msg.style.paddingTop= "0px";
    }, 3000);

};

function errorScript(err) {
    var err_msg = err;
    let h3_ele = document.querySelector('#err_msg')
    h3_ele.innerHTML=err_msg.replace(/"/g,"")
    h3_ele.style.color = "red";
    h3_ele.style.paddingTop = "10px";
    setTimeout(function() { 
        h3_ele.innerHTML=null;
        h3_ele.style.paddingTop = "0px";
    }, 10000); 
};

function showSpinner() {

    submitButton.style.visibility = "hidden";
    Spinner.style.visibility = "visible";

};

formSubmit.addEventListener('submit', function(e) {

    const allowed_ext = ['pdb', 'sdf'];

    if (Smiles.value) {
        showSpinner();
    } else if (allowed_ext.includes(MolFile.value.split('.').pop())) {
        showSpinner();
    } else {
        e.preventDefault();
        errorMessage();
    }

});

