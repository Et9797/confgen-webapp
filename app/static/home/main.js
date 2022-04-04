// Contact
$("#sendBtnContact").on("click", () => {
    $("#alertContact").css("visibility", "visible")
    setTimeout(() => $("#alertContact").css("visibility", "hidden"), 5000)
})

$("#contact-form").on("submit", (e) => {
    e.preventDefault()
    const formData = new FormData($("#contact-form")[0])
    fetch("/contact", {
        method: "POST",
        body: formData
    })
})

// mol ext can only be in separate files
$("#molRadio").on("click", () => {
    $("#separateFiles").prop("checked", true)
    $("#separateFiles").attr("onclick", 'return false;')
    $("#separateFiles").attr("onkeydown", 'return false;')
})

$("#pdbRadio, #sdfRadio").on("click", () => {
    $("#separateFiles").attr("onclick", '')
    $("#separateFiles").attr("onkeydown", '')
})

// On form submit
$("#main-form").on("submit", (e) => {
    e.preventDefault()
    const allowedExtensions = ["pdb", "sdf", "mol"]
    const formData = new FormData($("#main-form")[0])
    const smiles = Array.from(formData.entries())[1][1]
    const molFile = Array.from(formData.entries())[0][1]["name"]
    const noConfs = Array.from(formData.entries())[2][1]
   
    // check de vals en dan mainform.trigger(submit)

    const asyncSubmit = async () => {
        // Hide submit button, show generating button
        $(".submitBtn").css("visibility", "hidden")
        $(".generatingBtn").css("visibility", "visible")
        
        // Calls Flask /generate route to start task in the background
        const r = await fetch("/generate", {
            method: "POST",
            body: formData
        })
        const rJSON = await r.json()
        const uniq_id = rJSON["uniq_id"]
        const task_id = rJSON["task_id"]
        
        // Poll the status of the task every 5 seconds
        const timeout = (ms) => new Promise(resolve => setTimeout(resolve, ms))
        const polling = (async () => {
            while (true) {
                const r = await fetch(`/task_status/${task_id}`, {
                    method: "GET"
                })
                const status = await r.json()
                if (status["state"] == "SUCCESS") {
                    $(".generatingBtn").css("visibility", "hidden")
                    showAlert("success", noConfs, 10000)
                    $(".download-reset-btns").css("visibility", "visible")
                    $("#downloadBtn").prop("href", `/${uniq_id}`)
                    break
                } else if (status["state"] == "FAILURE") {
                    $(".generatingBtn").css("visibility", "hidden")
                    $(".submitBtn").css("visibility", "visible")
                    showAlert("danger", null, 7000)
                    break
                } 
                await timeout(1000)
            }
        })()
    }

    if (smiles) {
        asyncSubmit()
    }
    else if (allowedExtensions.includes(molFile.split(".").pop())) {
        asyncSubmit()
    }
    else {
        e.preventDefault()
        showAlert("danger", null, 7000)
    }
})

function showAlert(type, noConfs, duration) {
    const alertFade = $(`.alerts .alert-${type}`)
    alertFade.empty()
    if (type == "success") {
        alertFade.append(`<i class="fas fa-check"></i> Successfully generated ${noConfs} conformers.`)
    } else if (type == "danger") {
        alertFade.append(`<i class="fas fa-times"></i> Something went wrong in RDKit. Check if the provided file or SMILES is correct.`) 
    }
    alertFade.addClass("fadeInAlert")
    setTimeout(() => {
        alertFade.removeClass("fadeInAlert")
        alertFade.addClass("fadeOutAlert")
        setTimeout(() => {
            alertFade.removeClass("fadeOutAlert")
        }, 500)
    }, duration)
}