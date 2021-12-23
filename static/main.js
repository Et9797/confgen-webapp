// Contact
$(".sendBtnContactForm").on("click", () => {
    $(".alertContact").css("visibility", "visible")
    setTimeout(() => $(".alertContact").css("visibility", "hidden"), 5000)
})

$("#contact-form").on("submit", (e) => {
    e.preventDefault()
    const formData = new FormData($("#contact-form")[0])
    fetch("/contact", {
        method: "POST",
        body: formData
    })
})

// mol ext can only be in seperate files
$("#molRadio").on("click", () => {
    $("#seperateFiles").prop("checked", true)
    $("#seperateFiles").attr("onclick", 'return false;')
    $("#seperateFiles").attr("onkeydown", 'return false;')
})

$("#pdbRadio, #sdfRadio").on("click", () => {
    $("#seperateFiles").attr("onclick", '')
    $("#seperateFiles").attr("onkeydown", '')
})

// Do on form submit
$("#main-form").on("submit", (e) => {
    const allowedExtensions = ["pdb", "sdf", "mol"]
    const formData = new FormData($("#main-form")[0])
    const smiles = Array.from(formData.entries())[1][1]
    const molFile = Array.from(formData.entries())[0][1]["name"]
    const hideSubmitBtn = () => {
        $(".submitBtn").css("visibility", "hidden")
        $(".generatingBtn").css("visibility", "visible")
    }
    if (smiles) {
        hideSubmitBtn()
    }
    else if (allowedExtensions.includes(molFile.split(".").pop())) {
        hideSubmitBtn()
    }
    else {
        e.preventDefault()
        showAlert("danger", 5000, null, "Check if the provided file/SMILES is correct.")
    }
})

// Checks if download button is available on the DOM
if ($("#downloadBtn").length) {
    const noConfs = $("#downloadBtn").prop("href").split("=").slice(-1).toString()
    showAlert("success", 5000, noConfs, null)
}

// If error display alert
if ($(".Error").length) {
    const error = $(".Error").prop("id")
    switch (error) {
        case "rdkit":
            showAlert("danger", 10000, null, `An error occurred in RDKit when trying to generate conformers. Check if
                                              the provided file/SMILES is correct.`)
            break
        case "NIH":
            showAlert("danger", 10000, null, `An error occurred trying to convert PDB to Smiles, unable to fetch NIH Api`)
            break
    }
}

// Alerts
function showAlert(type, duration, noConfs, exception) {
    const alertFade = $(`.alerts .alert-${type}`)
    alertFade.empty()
    if (type == "success") {
        alertFade.append(`<i class="fas fa-check"></i> Successfully generated ${noConfs} conformers.`)
    } else if (type == "danger") {
        alertFade.append(`<i class="fas fa-times"></i> Something went wrong. ${exception}`) 
    }
    alertFade.css("visibility", "visible")
    alertFade.addClass("alert-message-slideIn")
    setTimeout(() => {
        alertFade.removeClass("alert-message-slideIn")
        alertFade.addClass("alert-message-slideOut")
        setTimeout(() => {
            alertFade.css("visibility", "hidden")
            alertFade.removeClass("alert-message-slideOut")
        }, 300)
    }, duration)
}
