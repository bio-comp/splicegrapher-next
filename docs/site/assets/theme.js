(() => {
  const STORAGE_KEY = "sgn-theme";

  function systemTheme() {
    return window.matchMedia("(prefers-color-scheme: dark)").matches ? "dark" : "light";
  }

  function getInitialTheme() {
    const stored = window.localStorage.getItem(STORAGE_KEY);
    if (stored === "light" || stored === "dark") {
      return stored;
    }
    return systemTheme();
  }

  function updateToggleLabels(theme) {
    const nextTheme = theme === "dark" ? "light" : "dark";
    const buttonText = nextTheme === "dark" ? "Dark mode" : "Light mode";
    document.querySelectorAll(".sgn-theme-toggle").forEach((button) => {
      button.setAttribute("aria-pressed", String(theme === "dark"));
      button.textContent = buttonText;
      button.setAttribute("title", `Switch to ${nextTheme} theme`);
    });
  }

  function applyTheme(theme) {
    document.documentElement.setAttribute("data-theme", theme);
    updateToggleLabels(theme);
  }

  function toggleTheme() {
    const current = document.documentElement.getAttribute("data-theme") || getInitialTheme();
    const next = current === "dark" ? "light" : "dark";
    window.localStorage.setItem(STORAGE_KEY, next);
    applyTheme(next);
  }

  document.addEventListener("DOMContentLoaded", () => {
    applyTheme(getInitialTheme());
    document.querySelectorAll(".sgn-theme-toggle").forEach((button) => {
      button.addEventListener("click", toggleTheme);
    });
  });
})();
