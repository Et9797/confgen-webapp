$(document).ready(() => {
    // Poll the status of the task every 2s
    const timeout = (ms) => new Promise(resolve => setTimeout(resolve, ms))
    const polling = (async () => {
        const params = new URLSearchParams(window.location.search)
        const TASK_ID = params.get("task_id")
        while (true) {
            const r = await fetch(`/task_status/${TASK_ID}`, {
                method: "GET"
            })
            const status = await r.json()
            if (status["state"] != "PENDING") {
                window.location.replace(`/results/${TASK_ID}`)
            }
            await timeout(2000)
        }
    })()
})