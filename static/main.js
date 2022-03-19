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

// On form submit
$("#main-form").on("submit", (e) => {
    // e.preventDefault()
    const allowedExtensions = ["pdb", "sdf", "mol"]
    const formData = new FormData($("#main-form")[0])
    const smiles = Array.from(formData.entries())[1][1]
    const molFile = Array.from(formData.entries())[0][1]["name"]
    const noConfs = Array.from(formData.entries())[2][1]
    const hideSubmitBtn = () => {
        // Hide submit button, show generating button
        $(".submitBtn").css("visibility", "hidden")
        $(".generatingBtn").css("visibility", "visible")
        
        // Call Flask /generate route to generate conformers in the background
        // fetch("/generate", {
        //     method: "POST",
        //     body: formData})
        // .then(response => response.json())
        // .then(jsonResponse => {
        //     // Download button
        //     $(".generatingBtn").css("visibility", "hidden")
        //     $(".submitBtn").css("visibility", "hidden")
        //     showAlert("success", 5000, noConfs, null)
        //     $(".download-reset-btns").css("visibility", "visible")
        //     $("#downloadBtn").prop("href", `/${jsonResponse["job_id"]}`)
        // })
        
        // .catch(exc => {
        //     showAlert("danger", 10000, null, exc.message)
        //     $(".submitBtn").css("visibility", "visible")
        //     $(".generatingBtn").css("visibility", "hidden")
        // })

    }

    if (smiles) {
        hideSubmitBtn()
    }
    else if (allowedExtensions.includes(molFile.split(".").pop())) {
        hideSubmitBtn()
    }
    else {
        e.preventDefault()
        showAlert("danger", 5000, null)
    }
})


// Polling status


// Checks if download button is available on the DOM
// if ($("#downloadBtn").length) {
//     const noConfs = $("#downloadBtn").prop("href").split("=").slice(-1).toString()
//     showAlert("success", 5000, noConfs, null)
// }



function showAlert(type, duration, noConfs) {
    const alertFade = $(`.alerts .alert-${type}`)
    alertFade.empty()
    if (type == "success") {
        alertFade.append(`<i class="fas fa-check"></i> Successfully generated ${noConfs} conformers.`)
    } else if (type == "danger") {
        alertFade.append(`<i class="fas fa-times"></i> Something went wrong in RDKit. Check if the provided file or SMILES is correct.`) 
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
