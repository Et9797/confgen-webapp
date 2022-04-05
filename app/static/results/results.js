$(document).ready(() => {
    // Poll the status of the task every 2s
    const timeout = (ms) => new Promise(resolve => setTimeout(resolve, ms))
    const polling = (async () => {
        const params = new URLSearchParams(window.location.search)
        const taskID = params.get("task_id")
        while (true) {
            const r = await fetch(`/task_status/${taskID}`, {
                method: "GET"
            })
            const status = await r.json()
            console.log(status);
            if (status["state"] == "SUCCESS") {
                $("#generating, #notification").css("display", "none")
                $("#finished").css("display", "flex")
                $(".addClass").toggleClass("animateElement")
                break
            } else if (status["state"] == "FAILURE") {
                $("#generating, #notification").css("display", "none")
                $("#error").css("display", "flex")
                $(".addClass").toggleClass("animateElement")
                $(".error-message").prop("innerText", status["info"])
                break
            }
            await timeout(2000)
        }
    })()
})