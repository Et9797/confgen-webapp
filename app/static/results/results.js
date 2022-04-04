// polling code
// uniq id en taskid

$(".addClass").toggleClass("animateElement");

//In this cases you should probably try .css('display','none') to hide and .css('display','') to show


// Poll the status of the task every 5 seconds
const timeout = (ms) => new Promise(resolve => setTimeout(resolve, ms))
const polling = (async () => {
    console.log("hi")
    // while (true) {
    //     await timeout(500)
    //     console.log("hello");
    // }
})()
    //     const r = await fetch(`/task_status/${task_id}`, {
    //         method: "GET"
    //     })
    //     const status = await r.json()
    //     console.log(status)
    //     console.log("hello");
    //     if (status["state"] == "SUCCESS") {
    //         $(".generatingBtn").css("visibility", "hidden")
    //         showAlert("success", noConfs, 10000)
    //         $(".download-reset-btns").css("visibility", "visible")
    //         $("#downloadBtn").prop("href", `/${uniq_id}`)
    //         break
    //     } else if (status["state"] == "FAILURE") {
    //         $(".generatingBtn").css("visibility", "hidden")
    //         $(".submitBtn").css("visibility", "visible")
    //         showAlert("danger", null, 7000)
    //         break
    //     } 
    //     await timeout(500)
    // }


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