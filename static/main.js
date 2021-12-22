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

$("#molRadio").on("click", () => {
    $("#seperateFiles").prop("checked", true)
    $("#seperateFiles").attr("onclick", 'return false;')
    $("#seperateFiles").attr("onkeydown", 'return false;')
})

$("#pdbRadio, #sdfRadio").on("click", () => {
    $("#seperateFiles").attr("onclick", '')
    $("#seperateFiles").attr("onkeydown", '')
})

$("#main-form").on("submit", (e) => {
    e.preventDefault()
    const allowedExtensions = ["pdb", "sdf", "mol"]
    const formData = new FormData($("#main-form")[0])
    const smiles = Array.from(formData.entries())[1][1]
    const molFile = Array.from(formData.entries())[0][1]["name"]
    const noConfs = Array.from(formData.entries())[2][1]
    const asyncSubmit = () => {
        $(".submitBtn").css("visibility", "hidden")
        $(".generatingBtn").css("visibility", "visible")
        fetch("/generate", {
            method: "POST",
            body: formData
        }).then(async (response) => {
            if (response.ok) {
                return response.json()
            }
            const error = await response.json()
            throw new Error(error.exception)
        }).then(jsonResponse => {
            $(".generatingBtn").css("visibility", "hidden")
            showAlert("success", 5000, noConfs, null)
            $(".download-reset-btns").css("visibility", "visible")
            $("#downloadBtn").prop("href", `/generate/${jsonResponse["job_id"]}`)
        }).catch(exc => {
            $(".generatingBtn").css("visibility", "hidden")
            $(".submitBtn").css("visibility", "visible")
            showAlert("danger", 10000, null, exc.message)
        })
    }
    if (smiles) {
        asyncSubmit()
    }
    else if (allowedExtensions.includes(molFile.split(".").pop())) {
        asyncSubmit()
    }
    else {
        showAlert("danger", 5000, null, "provide a valid file format.")
    }

})

function showAlert(type, duration, noConfs, exception) {
    const alertFade = $(`.alerts .alert-${type}`)
    alertFade.empty()
    if (type == "success") {
        alertFade.append(`<i class="fas fa-check"></i> Successfully generated ${noConfs} conformers.`)
    } else if (type == "danger") {
        alertFade.append(`<i class="fas fa-times"></i> Something went wrong. Check if the provided file/SMILES is correct. 
        Exception: ${exception}`) 
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
