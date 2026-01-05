/* ========================================
   Confgen Modern UI - JavaScript
   ======================================== */

// ============ Particle Animation ============
class ParticleNetwork {
  constructor(canvas) {
    this.canvas = canvas;
    this.ctx = canvas.getContext('2d');
    this.particles = [];
    this.particleCount = 60;
    this.maxDistance = 150;
    this.mousePosition = { x: null, y: null };

    this.resize();
    this.init();
    this.animate();

    window.addEventListener('resize', () => this.resize());
    window.addEventListener('mousemove', (e) => {
      this.mousePosition.x = e.clientX;
      this.mousePosition.y = e.clientY;
    });
  }

  resize() {
    this.canvas.width = window.innerWidth;
    this.canvas.height = window.innerHeight;
  }

  init() {
    this.particles = [];
    for (let i = 0; i < this.particleCount; i++) {
      this.particles.push({
        x: Math.random() * this.canvas.width,
        y: Math.random() * this.canvas.height,
        vx: (Math.random() - 0.5) * 0.5,
        vy: (Math.random() - 0.5) * 0.5,
        radius: Math.random() * 2 + 1,
        color: Math.random() > 0.5 ? '#06b6d4' : '#14b8a6'
      });
    }
  }

  animate() {
    this.ctx.clearRect(0, 0, this.canvas.width, this.canvas.height);

    // Update and draw particles
    this.particles.forEach((particle, i) => {
      // Move particle
      particle.x += particle.vx;
      particle.y += particle.vy;

      // Bounce off edges
      if (particle.x < 0 || particle.x > this.canvas.width) particle.vx *= -1;
      if (particle.y < 0 || particle.y > this.canvas.height) particle.vy *= -1;

      // Draw particle
      this.ctx.beginPath();
      this.ctx.arc(particle.x, particle.y, particle.radius, 0, Math.PI * 2);
      this.ctx.fillStyle = particle.color;
      this.ctx.globalAlpha = 0.6;
      this.ctx.fill();

      // Draw connections
      this.particles.slice(i + 1).forEach(other => {
        const dx = particle.x - other.x;
        const dy = particle.y - other.y;
        const distance = Math.sqrt(dx * dx + dy * dy);

        if (distance < this.maxDistance) {
          this.ctx.beginPath();
          this.ctx.moveTo(particle.x, particle.y);
          this.ctx.lineTo(other.x, other.y);
          this.ctx.strokeStyle = '#06b6d4';
          this.ctx.globalAlpha = 0.15 * (1 - distance / this.maxDistance);
          this.ctx.lineWidth = 1;
          this.ctx.stroke();
        }
      });

      // Mouse interaction
      if (this.mousePosition.x && this.mousePosition.y) {
        const dx = particle.x - this.mousePosition.x;
        const dy = particle.y - this.mousePosition.y;
        const distance = Math.sqrt(dx * dx + dy * dy);

        if (distance < 200) {
          this.ctx.beginPath();
          this.ctx.moveTo(particle.x, particle.y);
          this.ctx.lineTo(this.mousePosition.x, this.mousePosition.y);
          this.ctx.strokeStyle = '#22d3ee';
          this.ctx.globalAlpha = 0.2 * (1 - distance / 200);
          this.ctx.lineWidth = 1;
          this.ctx.stroke();
        }
      }
    });

    this.ctx.globalAlpha = 1;
    requestAnimationFrame(() => this.animate());
  }
}

// Initialize particle animation
document.addEventListener('DOMContentLoaded', () => {
  const canvas = document.getElementById('particleCanvas');
  if (canvas) {
    new ParticleNetwork(canvas);
  }
});

// ============ Modal Functions ============
function openModal(modalId) {
  const modal = document.getElementById(modalId);
  if (modal) {
    modal.classList.remove('hidden');
    modal.classList.add('active');
    // Prevent body scroll
    document.body.style.overflow = 'hidden';
  }
}

function closeModal(modalId) {
  const modal = document.getElementById(modalId);
  if (modal) {
    modal.classList.remove('active');
    setTimeout(() => {
      modal.classList.add('hidden');
    }, 300);
    document.body.style.overflow = '';
  }
}

// Close modal on overlay click
document.querySelectorAll('.modal-overlay').forEach(modal => {
  modal.addEventListener('click', (e) => {
    if (e.target === modal) {
      closeModal(modal.id);
    }
  });
});

// Close modal on Escape key
document.addEventListener('keydown', (e) => {
  if (e.key === 'Escape') {
    document.querySelectorAll('.modal-overlay.active').forEach(modal => {
      closeModal(modal.id);
    });
  }
});

// ============ Contact Form ============
$("#sendBtnContact").on("click", () => {
  $("#alertContact").removeClass("hidden");
  setTimeout(() => $("#alertContact").addClass("hidden"), 5000);
});

$("#contact-form").on("submit", (e) => {
  e.preventDefault();
  const formData = new FormData($("#contact-form")[0]);
  fetch("/contact", {
    method: "POST",
    body: formData
  });
});

// ============ File Upload / Drop Zone ============
const dropZone = document.getElementById('dropZone');
const fileInput = document.getElementById('fileInput');
const dropZoneContent = document.getElementById('dropZoneContent');
const fileSelected = document.getElementById('fileSelected');
const fileName = document.getElementById('fileName');
const removeFile = document.getElementById('removeFile');

if (dropZone) {
  // Click to upload
  dropZone.addEventListener('click', () => fileInput.click());

  // Drag events
  dropZone.addEventListener('dragover', (e) => {
    e.preventDefault();
    dropZone.classList.add('drag-over');
  });

  dropZone.addEventListener('dragleave', () => {
    dropZone.classList.remove('drag-over');
  });

  dropZone.addEventListener('drop', (e) => {
    e.preventDefault();
    dropZone.classList.remove('drag-over');

    const files = e.dataTransfer.files;
    if (files.length) {
      fileInput.files = files;
      showSelectedFile(files[0].name);
    }
  });

  // File input change
  fileInput.addEventListener('change', () => {
    if (fileInput.files.length) {
      showSelectedFile(fileInput.files[0].name);
    }
  });

  // Remove file
  removeFile.addEventListener('click', (e) => {
    e.stopPropagation();
    fileInput.value = '';
    hideSelectedFile();
  });
}

function showSelectedFile(name) {
  dropZoneContent.classList.add('hidden');
  fileSelected.classList.remove('hidden');
  fileName.textContent = name;
}

function hideSelectedFile() {
  dropZoneContent.classList.remove('hidden');
  fileSelected.classList.add('hidden');
}

// ============ Conformer Count Buttons ============
const decreaseBtn = document.getElementById('decreaseConf');
const increaseBtn = document.getElementById('increaseConf');
const confInput = document.getElementById('noConformers');

if (decreaseBtn && increaseBtn && confInput) {
  decreaseBtn.addEventListener('click', () => {
    const currentVal = parseInt(confInput.value) || 100;
    const newVal = Math.max(25, currentVal - 25);
    confInput.value = newVal;
  });

  increaseBtn.addEventListener('click', () => {
    const currentVal = parseInt(confInput.value) || 100;
    const newVal = Math.min(500, currentVal + 25);
    confInput.value = newVal;
  });
}

// ============ MOL Format Toggle Logic ============
$("#molRadio").on("click", () => {
  $("#separateFiles").prop("checked", true);
  $("#separateFiles").prop("disabled", true);
});

$("#pdbRadio, #sdfRadio").on("click", () => {
  $("#separateFiles").prop("disabled", false);
});

// ============ Form Validation & Submit ============
$("#main-form").on("submit", (e) => {
  const allowedExtensions = ["pdb", "sdf", "mol"];
  const formData = new FormData($("#main-form")[0]);
  const smiles = formData.get('SMILES');
  const molFile = formData.get('molFile');
  const molFileName = molFile ? molFile.name : '';

  // Hide previous error
  $("#errorFeedback").addClass("hidden");

  if (smiles && smiles.trim()) {
    // Show loading state
    $("#submitBtn").addClass("loading");
    return true;
  } else if (molFileName && allowedExtensions.includes(molFileName.split(".").pop().toLowerCase())) {
    // Show loading state
    $("#submitBtn").addClass("loading");
    return true;
  } else {
    e.preventDefault();
    $("#errorFeedback").removeClass("hidden");
    // Scroll to error
    $("#errorFeedback")[0].scrollIntoView({ behavior: 'smooth', block: 'center' });
    return false;
  }
});

// Make modal functions globally available
window.openModal = openModal;
window.closeModal = closeModal;
