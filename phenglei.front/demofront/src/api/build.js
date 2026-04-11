import axios from "axios";

const BASE_URL = "http://localhost:8080";

export const buildExe = () => axios.get(BASE_URL + "/build-exe");
export const runExe = () => axios.get(BASE_URL + "/run-exe");