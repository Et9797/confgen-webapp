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
    const allowedExtensions = ["pdb", "sdf", "mol"]
    const formData = new FormData($("#main-form")[0])
    const smiles = Array.from(formData.entries())[1][1]
    const molFile = Array.from(formData.entries())[0][1]["name"]
   
    if (smiles) {
        return true
    }
    else if (allowedExtensions.includes(molFile.split(".").pop())) {
        return true
    }
    else {
        e.preventDefault()
        $(".invalid-feedback").css("display", "flex").delay(5000).fadeOut("fast")
    }

})